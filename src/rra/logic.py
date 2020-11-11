from typing import Dict, List, Set, Generator, Tuple
from .data import DataProvider

from pandas import DataFrame


def add_standards(dp:DataProvider, compounds:List[int]) -> Set:
    compound_set = set(compounds)

    while True:
        standard_set = set()
        for coeq in [dp.eq_classes.get(dp.classifications.get(c)) for c in compound_set]:
            if coeq is None: continue
            standard_set = standard_set.union(set(coeq))
        standard_set = standard_set.difference(compound_set)
        if not standard_set:
            break
        compound_set = compound_set.union(standard_set)

    return compound_set


def get_compounds(df:DataFrame) -> Set[int]:
    return set(df['compound'].values)


def prediction_graph(
        dp:DataProvider,
        compound_set:Set[int]
    ) -> DataFrame:

    tag = 1
    predictions = dp.predictions_from(compound_set, tag=tag)
    additional_compounds = get_compounds(predictions).difference(compound_set)
    
    while additional_compounds:
        additional_compounds = add_standards(dp, additional_compounds)
    
        compound_set = compound_set.union(additional_compounds)
        
        tag += 1
        downstream_predictions = dp.predictions_from(additional_compounds, tag=tag)
        predictions = predictions.append(downstream_predictions, sort=False).drop_duplicates()
        
        additional_compounds = get_compounds(predictions).difference(compound_set)
        
    return predictions


def reaction_prediction(
        dp:DataProvider, reac_name:str
    ):

    reaction_ids = dp.reaction_id_from_name(reac_name)
    substrates, products = dp.compounds_of_reactions(reaction_ids)
    
    predictions = []
    for compounds in [substrates, products]:
        compounds_plus = add_standards(dp, get_compounds(compounds))
        predicted = prediction_graph(dp, compounds_plus)
        predictions.append(predicted)
    
    SUBSTRATES,PRODUCTS = 0,1
    return predictions[SUBSTRATES], predictions[PRODUCTS]


class LimitExceeded(Exception): pass

Reaction = int # reaction id
EqClass = int # eqcl id
Path = List[Reaction]
Match = List[Path]

def _path_finder(
        tree:DataFrame,
        starteducts:List[EqClass],
        endproduct:EqClass,
        stations:Dict[EqClass,List[Reaction]],
        hitting_limit:int,
        max_depth:int
    ) -> Generator[Path,None,None]:

    hitting = tree[
        (tree.eqcl == endproduct)
        & tree.isproduct
    ].reaction.drop_duplicates().values

    if hitting_limit is not None and len(hitting) > hitting_limit:
        raise LimitExceeded(f"hit limit of {hitting_limit} reached", endproduct, len(hitting)) 

    if not len(hitting): # the endproduct is not reachable
        return
    
    if max_depth is not None:
        max_depth -= 1
    
    for pred, reac in tree[
            tree.reaction.isin(hitting)
            & ~tree.isproduct
        ].apply(lambda row: (row.eqcl, row.reaction), axis=1)\
         .drop_duplicates()\
         .values:

        if pred in starteducts: # all done
            yield [reac]
            continue

        if max_depth is not None and max_depth < 0:
            continue
        
        if pred in stations:
            for path in stations[pred]:
                if not reac in path:
                    yield path + [reac]
            continue
        stations[pred] = []

        for path in _path_finder(tree, starteducts, pred,
                                 stations, hitting_limit, max_depth):
            stations[pred].append(path)
            if not reac in path:
                yield path + [reac]


def _reversed_path_finder(
        tree:DataFrame,
        starteduct:EqClass,
        endproducts:List[EqClass],
        stations:Dict[EqClass,List[Reaction]],
        hitting_limit:int,
        max_depth:int
    ) -> Generator[Path,None,None]:

    hitting = tree[
        (tree.eqcl == starteduct)
        & ~tree.isproduct
    ].reaction.drop_duplicates().values

    if hitting_limit is not None and len(hitting) > hitting_limit:
        raise LimitExceeded(f"hit limit of {hitting_limit} reached", starteduct, len(hitting)) 

    if not len(hitting): # the starteduct is inert
        return

    if max_depth is not None:
        max_depth -= 1

    for succ, reac in tree[
            tree.reaction.isin(hitting)
            & tree.isproduct
        ].apply(lambda row: (row.eqcl, row.reaction), axis=1)\
         .drop_duplicates()\
         .values:

        if succ in endproducts: # all done
            yield [reac]
            continue

        if max_depth is not None and max_depth < 0:
            continue
        
        if succ in stations:
            for path in stations[succ]:
                if not reac in path:
                    yield [reac] + path
            continue
        stations[succ] = []

        for path in _reversed_path_finder(tree, succ, endproducts,
                                          stations, hitting_limit, max_depth):
            stations[succ].append(path)
            if not reac in path:
                yield [reac] + path


def _filter_irrelevant(dp, df):
    return df[~df.eqcl.isin(dp.irrelevant_eq_classes)]

def _filter_transitions(dp, substrates, products):
    subs, pros = substrates[:], products[:]
    for tr1, tr2 in dp.transitions:
        if tr1 in subs and tr2 in pros:
            subs.remove(tr1)
            pros.remove(tr2)
        if tr2 in subs and tr1 in pros:
            subs.remove(tr2)
            pros.remove(tr1)
    return subs, pros

def get_match_compounds(dp:DataProvider, reac_ids:List[Reaction]) -> Tuple[List[EqClass], List[EqClass]]:
    subdf, prodf = dp.compounds_of_reactions(reac_ids)
    subdf, prodf = map(lambda df: _filter_irrelevant(dp, df), (subdf, prodf))
    return _filter_transitions(dp, list(subdf.eqcl.values), list(prodf.eqcl.values))


_FORWARD_ = True
_REVERSE_ = False
def pathes_of_reaction(dp:DataProvider, 
        sub_eqcls:List[EqClass],
        prd_eqcls:List[EqClass],
        forward_hitlimit=20,
        backward_hitlimit=20,
        step_max=None,
        suppress_warnings=False
    ) -> Generator[Path,None,None]:

    hitlimit_exceeding = {}
    for prd_eqcl in prd_eqcls:
        try:
            for path in _path_finder(dp.fulltree, sub_eqcls, prd_eqcl, dict(), hitting_limit=forward_hitlimit, 
            max_depth=step_max):
                # even though the step_max was passed over to _path_finder, it still needs to be used as
                # a filter, since the latter does not truncate pathes by itself.
                if step_max is None or step_max >= len(path):
                    yield _FORWARD_, path
        except LimitExceeded as le:
            hitlimit_exceeding[le.args[1]] = le.args[2]
    if hitlimit_exceeding:
        if not suppress_warnings:
            from sys import stderr
            stderr.write(f"some of the products have a very high (>{forward_hitlimit}) hit rate: {hitlimit_exceeding}\n")
            stderr.write("they were skipped in the search. the search is repeated from the other end.\n")
            stderr.write("this may lead to duplicated pathes.\n")
        for sub_eqcl in sub_eqcls:
            for path in _reversed_path_finder(dp.fulltree, sub_eqcl, prd_eqcls, dict(), hitting_limit=backward_hitlimit,
            max_depth=step_max):
                if step_max is None or step_max >= len(path):
                    yield _FORWARD_, path

    hitlimit_exceeding = {}
    for sub_eqcl in sub_eqcls:
        try:
            for path in _path_finder(dp.fulltree, prd_eqcls, sub_eqcl, dict(), hitting_limit=forward_hitlimit,
            max_depth=step_max):
                if step_max is None or step_max >= len(path):
                    yield _REVERSE_, path
        except LimitExceeded as le:
            hitlimit_exceeding[le.args[1]] = le.args[2]
    if hitlimit_exceeding:
        if not suppress_warnings:
            from sys import stderr
            stderr.write(f"some of the substrates have a very high (>{forward_hitlimit}) hit rate: {hitlimit_exceeding}\n")
            stderr.write("they were skipped in the search. the search is repeated from the other end\n")
            stderr.write("this may lead to duplicated pathes\n")
        for prd_eqcl in prd_eqcls:
            for path in _reversed_path_finder(dp.fulltree, prd_eqcl, sub_eqcls, dict(), hitting_limit=backward_hitlimit,
            max_depth=step_max):
                if step_max is None or step_max >= len(path):
                    yield _REVERSE_, path


def match_reaction(
        dp:DataProvider,
        reac_id:Reaction,
        forward_hitlimit:int=20,
        backward_hitlimit:int=20,
        step_max:int=None,
        supermatch_allowed:bool=True,
        suppress_warnings=False
) -> Generator[Match, None, None]:

    substrates, products = get_match_compounds(dp, [reac_id])

    allpathes = []
    for sense, path in pathes_of_reaction(dp, substrates, products, forward_hitlimit, backward_hitlimit, step_max, suppress_warnings):
        if (sense, path) not in allpathes: # here python compares the values and not the reference
            allpathes.append((sense, path))
    # allpathes is now a non-redundant list of pathes
    # leading from at least one substrate to at least one product of the reaction.

    original_reaction = set(substrates + products) if supermatch_allowed else None
    for match in _match_finder(allpathes, substrates, products, [], dp, original_reaction):
        yield match


_END_ = True
_START_ = False
def _match_finder(allpathes, substrates, products, match_attempt, dp, original_reaction):
    while True:
        try:
            sense, frst_path = allpathes.pop()
        except IndexError: # list is empty
            return
        
        subs = substrates[:]
        prods = products[:]
        
        # check for common path tails that split off at some point and never come back
        endintegration = _integrates_at_the_end(dp, match_attempt, sense, frst_path)
        startintegration = _integrates_from_the_start(dp, match_attempt, sense, frst_path)

        if endintegration and not startintegration:
            fine = _pop_substances(dp, prods, sense, frst_path, _END_, original_reaction)
            # fine: truly integrates towards the end, i.e. the end products of frst_path were still listed (forward)
        if startintegration and not endintegration:
            fine = _pop_substances(dp, subs, sense, frst_path, _START_, original_reaction)
            # fine: truly integrates from the start, i.e. the start substrates of frst_path were still listed (forward)
        if not endintegration and not startintegration:
            fine = _pop_substances(dp, prods, sense, frst_path, _END_, original_reaction) \
               and _pop_substances(dp, subs, sense, frst_path, _START_, original_reaction)
            # fine: products _and_ substrates were still listed
        if endintegration and startintegration:
            fine = _pop_substances(dp, prods, sense, frst_path, _END_, original_reaction) \
                or _pop_substances(dp, subs, sense, frst_path, _START_, original_reaction)
            # fine: at least one side, products or substrates, was still listed
        
        if not fine: # cannot reduce substrate nor product number
            continue
        
        next_match_attempt = match_attempt[:]
        next_match_attempt.append((sense, frst_path))

        if not subs and not prods: # all reaction partners have been assigned to a path now
            yield next_match_attempt
        
        else: # not all reaction partners have been assigned to a path yet
            for match in _match_finder(allpathes[:], subs, prods, next_match_attempt, dp, original_reaction):
                yield match


def _pop_substances(dp, slist, sense, path, side, original_reaction):
    # NOTE! (here as well as anywhere)
    # stoichiometry is not taken into account nowhere in the rra.logic module
    # if the reaction has more than 1 substance per equivalent class
    # all of them are taken out at once if the path is ending there
    if (sense == _FORWARD_ and side == _END_) \
    or (sense == _REVERSE_ and side == _START_):
        prodf = dp.products_of_reaction(path[-1])
        substances = prodf.eqcl.drop_duplicates().values
    else:
        subdf = dp.substrates_of_reaction(path[0])
        substances = subdf.eqcl.drop_duplicates().values
    for s in substances:
        if s not in slist:
            if original_reaction is None or s in original_reaction:
                return False
    for s in substances:
        if original_reaction is None or s in original_reaction:
            slist.remove(s)
    return True


def _integrates_at_the_end(dp, match_attempt, sense, path):
    for attsense, attpath in match_attempt:
        # the direction of the pathes must match
        if attsense != sense:
            continue
        
        # find common start sequence
        i = 0
        while i < len(path) and i < len(attpath):
            if path[i] != attpath[i]:
                break
            i += 1
        # no common start sequence - no integration
        if i == 0:
            continue
            
        # The last common reaction, path[i-1] (=attpath[i-1]), must be splitting 
        # and its products either substrates of path[i] or attpath[i], not both.
        # => intersection of substrates from path[i] and attpath[i] and products of path[i-1] must be empty.
        subdf, prodf = dp.compounds_of_reaction(path[i-1])
        products = set(prodf.eqcl.values) if sense == _FORWARD_ else set(subdf.eqcl.values)
        if i < len(path):
            subdf, prodf = dp.compounds_of_reaction(path[i])
            pathsubstrates = set(subdf.eqcl.values) if sense == _FORWARD_ else set(prodf.eqcl.values)
        else:
            pathsubstrates = set()
        if i < len(attpath):
            subdf, prodf = dp.compounds_of_reaction(attpath[i])
            attpathsubstrates = set(subdf.eqcl.values) if sense == _FORWARD_ else set(prodf.eqcl.values)
        else:
            attpathsubstrates = set()
        
        if (i >= len(path) and len(products) >= len(attpathsubstrates)) \
        or (i >= len(attpath) and len(products) >= len(pathsubstrates)) \
        or products.intersection(pathsubstrates).intersection(attpathsubstrates):
            continue

        # check whether the pathes merge again after the split
        mergeagain = False
        for step in path[i:]:
            if step in attpath[i:]:
                mergeagain = True
                break
        # if they do - no integration
        if mergeagain is True:
            continue
        
        # all tests passed
        return True

    return False


def _integrates_from_the_start(dp, match_attempt, sense, path):
    for attsense, attpath in match_attempt:
        # the direction of the pathes must match
        if attsense != sense:
            continue

        # find common end sequence
        i = -1
        while 1-i < len(path) and i < len(attpath):
            if path[i] != attpath[i]:
                break
            i -= 1
        # no common end sequence - no integration
        if i == -1:
            continue

        # The first common reaction, path[i+1] (=attpath[i+1]), must be splitting 
        # and its substrates either products of path[i] or attpath[i], not both.
        # => intersection of products from path[i] and attpath[i] and substrates of path[i+1] must be empty.
        subdf, prodf = dp.compounds_of_reaction(path[i+1])
        substrates = set(subdf.eqcl.values) if sense == _FORWARD_ else set(prodf.eqcl.values)
        if -i <= len(path):
            subdf, prodf = dp.compounds_of_reaction(path[i])
            pathproducts = set(prodf.eqcl.values) if sense == _FORWARD_ else set(subdf.eqcl.values)
        else:
            pathproducts = set()
        if -i <= len(attpath):
            subdf, prodf = dp.compounds_of_reaction(attpath[i])
            attpathproducts = set(prodf.eqcl.values) if sense == _FORWARD_ else set(subdf.eqcl.values)
        else:
            attpathproducts = set()
        # no true split - no integration
        if (-i > len(path) and len(substrates) >= len(attpathproducts)) \
        or (-i > len(attpath) and len(substrates) >= len(pathproducts)) \
        or substrates.intersection(pathproducts).intersection(attpathproducts):
            continue

        # check whether the pathes merge again after the split
        mergeagain = False
        for step in path[:i]:
            if step in attpath[:i]:
                mergeagain = True
                break
        # if the don't it's a true integration
        if mergeagain is True:
            continue
        
        # all tests passed
        return True

    return False


def sufficient_rule_sets(dp, matches, coln='btrule'):
    # collect unique single steps that combine to a full match
    stepsets = set()
    for match in matches:
        stepset = set()
        direction = None
        for path in match:
            direction = 'f' if path[0] and direction != 'r' \
                    else 'r' if not path[0] and direction != 'f' \
                    else 'b'
            stepset = stepset.union(set(path[1]))
        if direction:
            stepsets.add(direction+",".join([str(i) for i in sorted(list(stepset))]))
        else:
            raise Exception("unexpected empty match (unnecessary raising though)")

    # collect combinations of rules that trigger these steps
    def lmerge(a,b):
        for x in a:
            for y in set(b):
                yield x + [y]

    rulesets = set()
    for stepset in stepsets:
        ruleset = [[]]
        steps = [int(s) for s in stepset[1:].split(',')]
        for step in steps:
            rules = dp.rule_reaction.loc[dp.rule_reaction.reaction == step, coln]
            ruleset = list(lmerge(ruleset, rules))
        for rulcomb in ruleset:
            rulesets.add(stepset[0]+",".join(sorted(list(set([str(x) for x in rulcomb])))))

    return [(r[0],r[1:].split(',')) for r in rulesets]


def non_reducible_rule_sets(ample_rule_sets):
    nonreducibles = {'f':[], 'r':[], 'b':[]}
    for (direction,rule_set) in sorted(ample_rule_sets, key=lambda l:len(l[1])):
        if not any([all([nre in rule_set for nre in nr]) for nr in nonreducibles[direction]]):
            nonreducibles[direction].append(rule_set)
    return [('f', rule_set) for rule_set in nonreducibles['f']] \
         + [('r', rule_set) for rule_set in nonreducibles['r']] \
         + [('b', rule_set) for rule_set in nonreducibles['b']]

