
from typing import List, Set, Tuple, Dict

from pandas import DataFrame
from pandas import read_sql
from sqlalchemy import create_engine
from sqlalchemy.engine.base import Engine


def eq_classify(
        standards:DataFrame
    ) -> Tuple[Dict[int,Set[int]], Dict[int,int]]:

    def classify(eqc:Dict[int,Set[int]], clf:Dict[int,int], a:int, b:int):
        ac = clf.get(a, a)
        bc = clf.get(b, b)
        if ac in eqc and bc in eqc:
            eqc[ac] = eqc[ac].union(eqc.pop(bc))
            clf[bc] = ac
        elif ac in eqc:
            eqc[ac].add(bc)
            clf[bc] = ac
        elif bc in eqc:
            eqc[bc].add(ac)
            clf[ac] = bc
        else:
            eqc[ac] = set([a,b])
            clf[b] = ac
            clf[a] = ac

    eq_classes = {}
    classifications = {}
        
    standards.apply(lambda row: classify(eq_classes, classifications, row['compound'], row.standard), axis=1)
    
    return eq_classes, classifications


class DataProvider():
    def conjure(self, driver:str, host:str, dbname:str, user:str, magic:str) -> Engine:
        try:
            pws = [x.strip() for x in open(magic).readlines()]
        except FileNotFoundError:
            import os
            raise Exception(f"no {magic} in {os.getcwd()}")
        return create_engine(f'{driver}://{user}:{pws[1]}@{host}/{dbname}')

    def enhance_rrole(self, rrole:DataFrame, standards:DataFrame):
        self.eq_classes, self.classifications = eq_classify(standards)
        rrole['eqcl'] = rrole.apply(
            lambda row: self.classifications.get(row['compound'], row['compound']),
            axis=1
        )
        return rrole
    
    def remove_automorphisms(self, rrole:DataFrame):
        mcounter = rrole.groupby(['reaction']).count()
        monos = mcounter[mcounter.eqcl==2].index
        acounter = rrole[rrole.reaction.isin(monos)]\
                .loc[:,['reaction', 'eqcl', 'isproduct']]\
                .drop_duplicates()\
                .groupby(['reaction', 'eqcl'])\
                .count()
        autos = [i[0] for i in acounter[acounter.isproduct==2].index]
        return rrole[~rrole.reaction.isin(autos)], rrole[rrole.reaction.isin(autos)]

    def __init__(self, magic=None,
            reactor_id:int=None,
            rule_ids:List[int]=None,
            standardizer_ids:List[int]=None,
            reload=None,
            transition_file=None,
            driver='postgresql+psycopg2',
            host='imsb-ra-monitor.ethz.ch',
            dbname='reactions',
            user='pps',
            password=None
        ):
        
        if password:
            self.rradb = create_engine(f'{driver}://{user}:{password}@{host}/{dbname}')
        elif magic:
            self.rradb = self.conjure(driver, host, dbname, user, magic)
        elif not reload:
            raise Exception('must be one of password, magic, reload')

        if reload:
            if any([rule_ids, reactor_id, standardizer_ids]):
                from sys import stderr
                stderr.write('reloading data from files, rule_ids, reactor_id and standardizer_ids paramter will be ignored')
            self.reload(reload)
        
        else:
            # cache data
            std_constraint = f" and standardizer_id in ({','.join(standardizer_ids)})" if standardizer_ids else ''
            self.standards = read_sql(f"select * from standards where standard != compound{std_constraint}", self.rradb)

            self.rrole = self.enhance_rrole(
                read_sql("select * from rrole", self.rradb),
                self.standards
            )
            self.rrole, self.srole = self.remove_automorphisms(self.rrole)
            self.rrole.set_index('reaction', inplace=True)

            rr_constraint = f" and rule in ({','.join([str(rid) for rid in rule_ids])}) and reactor = {reactor_id}" if (rule_ids and reactor_id)\
                    else f" and rule in ({','.join([str(rid) for rid in rule_ids])})" if rule_ids \
                    else f" and reactor = {reactor_id}" if reactor_id \
                    else ''
            self.rule_reaction = read_sql(f"""
                select reaction, rule, reactor, container btrule, simple simpfix
                from rule_reaction join rule on rule.id = rule
                where not excluded is TRUE{rr_constraint}""", self.rradb)
            self.fulltree = self.rule_reaction.loc[:,['reaction', 'rule', 'reactor']].merge(self.rrole, on='reaction')

            self.irrelevant_eq_classes = list(map(
                lambda c: self.classifications.get(c, c),
                read_sql("select id as cmpd from irrelevant", self.rradb).cmpd.values
            ))
            
            self.transitions = self.read_transitions_from_file_or_db(transition_file)
        
    def read_transitions_from_file_or_db(self, transition_file):
        trs = []
        from pandas import read_csv

        if transition_file:
            trs_df = read_csv(transition_file, header=None, index_col=None, sep="\t")
            trs_df.columns = ['a', 'b']
        else:
            trs_df = read_sql("select a, b from transition", self.rradb)
        
        def read_transition(x,y):
            try:
                (nx, ny) =  (int(x), int(y))
            except ValueError:
                try:
                    (nx, ny) = map(
                        lambda c: read_sql(
                            f"select id as cmpd from nsmiles where smiles = '{c}'", self.rradb
                        ).cmpd.values[0],
                        (x,y)
                    )
                except Exception as e:
                    raise e
            return self.classifications.get(nx, nx), self.classifications.get(ny, ny)
        trs_df.apply(lambda row: trs.append(read_transition(row.a, row.b)), axis=1)
        return trs


    def reaction_id_from_name(self, name:str) -> List[int]:
        return read_sql(f'''
            select reaction
            from edb_reaction
            where id = '{name}'
            ''', self.rradb).reaction.values

    def compounds_of_reactions(self, reaction_ids:List[int]) -> DataFrame:
        rrdf = self.rrole.loc[reaction_ids,:]
        return (
            rrdf[~rrdf.isproduct],
            rrdf[rrdf.isproduct]
        )

    def products_of_reaction(self, reaction_id:int) -> DataFrame:
        rrdf = self.rrole.loc[reaction_id,:]
        return rrdf[rrdf.isproduct]

    def substrates_of_reaction(self, reaction_id:int) -> DataFrame:
        rrdf = self.rrole.loc[reaction_id,:]
        return rrdf[~rrdf.isproduct]

    def compounds_of_reaction(self, reaction_id:int) -> Tuple[DataFrame]:
        rrdf = self.rrole.loc[reaction_id,:]
        return (
            rrdf[~rrdf.isproduct],
            rrdf[rrdf.isproduct]
        )

    def standards_of_compounds(self, compounds:Set[int]) -> DataFrame:
        down = self.standards[self.standards['compound'].isin(compounds)]
        down.columns = ['compound', 'standardizer', 'std_link']
        down = down.copy()
        down['dir'] = 'down'
        
        down = down.copy()
        up = self.standards[self.standards['standard'].isin(compounds)].loc[:,['standard','standardizer','compound']]
        up.columns = ['compound', 'standardizer', 'std_link']
        up = up.copy()
        up['dir'] = 'up'
    
        return down.append(up, sort=False)
    
    def predictions_from(self,
            compound_ids:Set[int],
            tag:int=None
        ) -> DataFrame:
        reaction_df = self.rrole[self.rrole['compound'].isin(compound_ids) & (~self.rrole.isproduct)]
        reactions = reaction_df.index.values
        df = self.rule_reaction[self.rule_reaction.reaction.isin(reactions)] \
            .merge(self.rrole.loc[reactions,:], on='reaction')
        if tag: df['tag'] = tag
        return df

    def dump(self, directory:str, overwrite_existing=True):
        if overwrite_existing:
            from os import mkdir
            try:
                mkdir(directory)
            except FileExistsError:
                pass

        self.rule_reaction.to_csv(f'{directory}/rule_reaction.csv', sep="\t", index=None)
        self.rrole.to_csv(f'{directory}/rrole.csv', sep="\t")
        self.srole.to_csv(f'{directory}/srole.csv', sep="\t", index=None)
        self.standards.to_csv(f'{directory}/standards.csv', sep="\t", index=None)
        
        with open(f'{directory}/eq_classes.txt', 'w') as deq:
            deq.write(f'{self.eq_classes}')
        with open(f'{directory}/classifications.txt', 'w') as deq:
            deq.write(f'{self.classifications}')
        with open(f'{directory}/irrelevant_cmpds.txt', 'w') as deq:
            deq.write(f'{self.irrelevant_eq_classes}')
        with open(f'{directory}/transitions.txt', 'w') as deq:
            deq.write(f'{self.transitions}')
    
    def reload(self, directory:str):
        from pandas import read_csv
        self.rule_reaction = read_csv(f'{directory}/rule_reaction.csv', sep="\t")
        self.rrole = read_csv(f'{directory}/rrole.csv', sep="\t").set_index('reaction')
        self.srole = read_csv(f'{directory}/srole.csv', sep="\t")
        self.standards = read_csv(f'{directory}/standards.csv', sep="\t")
        
        from ast import literal_eval
        self.eq_classes = literal_eval(open(f'{directory}/eq_classes.txt').read())
        self.classifications = literal_eval(open(f'{directory}/classifications.txt').read())
        self.irrelevant_eq_classes = literal_eval(open(f'{directory}/irrelevant_cmpds.txt').read())
        self.transitions = literal_eval(open(f'{directory}/transitions.txt').read())

        self.fulltree = self.rule_reaction.merge(self.rrole, on='reaction')

