
'''
Created on May 19, 2016

@author: me
'''
from re import split
from sqlalchemy import MetaData, Table, Column, UniqueConstraint, Integer, String, Boolean, ForeignKey, create_engine
from sqlalchemy.exc import ResourceClosedError, StatementError, IntegrityError
from sqlalchemy.orm.session import sessionmaker
Session = sessionmaker()

class DataBase(object):
    '''
    takes care of the data base
    '''
    engine = None
    metadata = None
    hand = None
    player = None
    action = None

    def __init__(self, params=None):
        '''
        creates and initializes the data base
        '''
        self.metadata = MetaData()
        self.engine = self._engine(params if params else {})
        self.nsmiles = Table('nsmiles', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('smiles', String(512), unique=True))
        self.alsmiles = Table('alsmiles', self.metadata,
                          Column('id', Integer(), ForeignKey('nsmiles.id')),
                          Column('smiles', String(512), index=True))
        self.progress = Table('progress', self.metadata,
                          Column('id', Integer(), ForeignKey('nsmiles.id'), unique=True),
                          Column('standardized', Boolean()),
                          Column('propagated', Boolean()))
        self.irrelevant = Table('irrelevant', self.metadata,
                          Column('id', Integer(), ForeignKey('nsmiles.id'), unique=True))
        self.transition = Table('transition', self.metadata,
                          Column('a', Integer(), ForeignKey('nsmiles.id'), unique=True),
                          Column('b', Integer(), ForeignKey('nsmiles.id'), unique=True))
        self.reaction = Table('reaction', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('signature', String(64), unique=True))
        self.rrole = Table('rrole', self.metadata,
                          Column('compound', Integer(), ForeignKey('nsmiles.id')),
                          Column('reaction', Integer(), ForeignKey('reaction.id')),
                          Column('isproduct', Boolean()))
        self.standardizer = Table('standardizer', self.metadata,
                          Column('id', Integer(), primary_key=True),
                          Column('name', String(32), unique=True, nullable=False),
                          Column('description', String(1024)))
        self.standards = Table('standards', self.metadata,
                          Column('compound', Integer(), ForeignKey('nsmiles.id')),
                          Column('standardizer', Integer(), ForeignKey('standardizer.id')),
                          Column('standard', Integer(), ForeignKey('nsmiles.id')),
                          UniqueConstraint('compound', 'standardizer'))
        self.rule = Table('rule', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('container', String(16), index=True),
                          Column('simple', String(16), unique=True),
                          Column('description', String(1024)))
        self.rule_exclusion = Table('rule_exclusion', self.metadata,
                          Column('rule', Integer(), ForeignKey('rule.id'), nullable=False, unique=True),
                          Column('reason', String(16)))
        self.edb_reaction = Table('edb_reaction', self.metadata,
                          Column('id', String(8), primary_key=True),
                          Column('reaction', Integer(), ForeignKey('reaction.id')),
                          Column('envipath_url', String(48)))
        self.edb_reaction_exclusion = Table('edb_reaction_exclusion', self.metadata,
                          Column('edb_reaction', String(8), ForeignKey('edb_reaction.id'), nullable=False, unique=True),
                          Column('reason', String(8)))
        self.edb_rule_reaction = Table('edb_rule_reaction', self.metadata,
                          Column('edb_reaction', String(8), ForeignKey('edb_reaction.id')),
                          Column('bt_rule', String(8)),
                          Column('step_no', Integer()))
        self.reactor = Table('reactor', self.metadata,
                          Column('id', Integer(), primary_key=True),
                          Column('name', String(16), unique=True))
        self.rule_reaction = Table('rule_reaction', self.metadata,
                          Column('reaction', Integer(), ForeignKey('reaction.id')),
                          Column('rule', Integer(), ForeignKey('rule.id')),
                          Column('reactor', Integer(), ForeignKey('reactor.id')),
                          Column('excluded', Boolean()))
        self.match_es = Table('match_es', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('reactor', Integer(), ForeignKey('reactor.id')),
                          Column('reaction', Integer(), ForeignKey('reaction.id')))
        self.match_parts = Table('match_parts', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('match', Integer(), ForeignKey('match_es.id')),
                          Column('final_compound', Integer(), ForeignKey('nsmiles.id')),
                          Column('initial_compound', Integer(), ForeignKey('nsmiles.id')),
                          Column('initial_rule', Integer(), ForeignKey('rule.id')),
                          Column('steps', Integer()),
                          Column('standardized', Boolean()))
        self.match_rules = Table('match_rules', self.metadata,
                          Column('match_part', Integer(), ForeignKey('match_parts.id')),
                          Column('rule', Integer(), ForeignKey('rule.id')),
                          Column('step', Integer()))
        self.match_standardizers = Table('match_standardizers', self.metadata,
                          Column('match_part', Integer(), ForeignKey('match_parts.id')),
                          Column('standardizer', Integer(), ForeignKey('standardizer.id')))
        self.generationsize = Table('generationsize', self.metadata,
                          Column('id', Integer(), primary_key=True, autoincrement=True),
                          Column('reactor', Integer(), ForeignKey('reactor.id')),
                          Column('reaction', Integer(), ForeignKey('reaction.id')),
                          Column('generation', Integer()),
                          Column('nproducts', Integer()))
        self.metadata.create_all(self.engine)
        self.connection = self.engine.connect()
        self.sess = Session(bind=self.connection)

    def insert(self, row, reactorid=None):
        '''
        expects a line in this format:
        <SMILES[.SMILES]>\\t<label|\\[label[, \[label\]]>\\t<SMILES[.SMILES]>
        '''
        s,l,p = row.strip().split("\t")
        substrates = s.split(".")
        products = p.split(".")
        labels = split(r'\s*,\s*',l.strip(r"[\[\]]")) if l[0]=="[" and l[-1]=="]" else [l]

        for compound in substrates + products:
            try:
                self._execute(self.nsmiles.insert(), {'smiles':compound})
            except IntegrityError:
                # at first this wasn't necessary here only later when inserting reactions
                # perhaps because of the order session and connection are used to interact
                # with the database in a mixed mode?
                # TODO: stick to session and get rid of the connection in self._execute
                self.sess.commit()

        subids = [_[0] for _ in self.sess.query(self.nsmiles.c.id).filter(self.nsmiles.c.smiles.in_(substrates)).all()]
        prodids = [_[0] for _ in self.sess.query(self.nsmiles.c.id).filter(self.nsmiles.c.smiles.in_(products)).all()]
        signature = self.signature(subids, prodids)
        #print (signature)

        skiproles = False
        try:
            self._execute(self.reaction.insert(), {'signature':signature})
        except IntegrityError:
            # seems to be necessary (or at least a hack) to avoid an InternalError
            self.sess.commit()
            skiproles = True
        reacid = self.sess.query(self.reaction.c.id).filter(self.reaction.c.signature == signature).one()[0]
        #print (reacid)

        if not skiproles:
            for subid in subids:
                self._execute(self.rrole.insert(), {'compound': subid, 'reaction': reacid, 'isproduct': False})
            for prodid in prodids:
                self._execute(self.rrole.insert(), {'compound': prodid, 'reaction': reacid, 'isproduct': True})

        for label in labels:
            if reactorid:
                ruleid = self.extract_rule(label)
                if not self.sess.query(self.rule_reaction).\
                        filter(self.rule_reaction.c.reactor == reactorid).\
                        filter(self.rule_reaction.c.rule == ruleid).\
                        filter(self.rule_reaction.c.reaction == reacid).count():
                    self._execute(self.rule_reaction.insert(), {'reactor':reactorid, 'rule':ruleid, 'reaction':reacid})
            else:
                if not self.sess.query(self.edb_reaction).\
                        filter(self.edb_reaction.c.envipath_url == label).count():
                    autokey = "g{0:04d}".format(
                            self.sess.query(self.edb_reaction).filter(self.edb_reaction.c.id.like('g%')).count() + 1
                    )
                    self._execute(self.edb_reaction.insert(), {'reaction':reacid,'envipath_url':label,'id':autokey})
                else:
                    self._execute(self.edb_reaction.update().\
                            where(self.edb_reaction.c.envipath_url == label).\
                            values(reaction = reacid))

        self.sess.commit()

    def _execute(self, *args):
        try:
            return self.connection.execute(*args)
        except StatementError as se:
            if se.orig.__class__ is ResourceClosedError:
                self.connection = self.engine.connect()
                return self._execute(*args)
            else:
                raise

    def _engine(self, params):
        dbtype = params.get('dbtype', 'postgres')
        if dbtype == 'sqlite':
            enginestr = "sqlite:///%s" % params.get('dbname', ':memory:')
        elif dbtype == 'postgres':
            enginestr = "postgresql+psycopg2://%s%s@%s:%d/%s" % (
                params.get('username', 'poker'),
                (":%s" % params.get('password')) if 'password' in params.keys() else '',
                params.get('host', 'localhost'),
                params.get('port', 5432),
                params.get('dbname', 'poker'))
        else:
            raise ValueError("not yet able to connect to %s databases" % dbtype)
        return create_engine(enginestr)

    def query(self, statement):
        return self._execute(statement).fetchall()

    def signature(self, subids, prodids):
        subids.sort()
        prodids.sort()
        return ".".join([str(_) for _ in subids])+">>"+".".join([str(_) for _ in prodids])

    def extract_rule(self, label):
        try: simple = int(label)
        except:
            simple = label.split("-")[-1]
        try:
            return self.sess.query(self.rule.c.id).filter(self.rule.c.simple == str(simple)).one()[0]
        except:
            print("extraction failure for simple rule %s" % simple)
            raise

    def import_rule(self, row):
        '''
        expects a line in this format:
        <btid>\t<rxnid>
        '''
        btid,rxnid = row.strip().split(",")
        try:
            self._execute(self.rule.insert(), {'simple':rxnid, 'container':btid})
        except:
            print("failed to insert ", btid, rxnid)
        #self.sess.commit()    

    def map_edb_rule_reaction(self, reacid, btid, rank):
        try:
            self._execute(self.edb_rule_reaction.insert(),
                {'edb_reaction': reacid,
                 'bt_rule': btid,
                 'step_no': rank})
        except IntegrityError:
            print("failed to insert ", reacid, btid, rank)
        except:
            print("failed to insert ", reacid, btid, rank)
            raise
        self.sess.commit()    

    def map_edb_reaction(self, reacid, mapid):
        try:
            self._execute(self.edb_reaction.insert(),
                {'envipath_url': reacid,
                 'id': mapid})
        except IntegrityError:
            print("failed to insert ", reacid, mapid)
        except:
            print("failed to insert ", reacid, mapid)
            raise
        self.sess.commit()

    def insert_standard(self, standardizerid, row):
        '''
        expects a line in this format:
        <int>\t<SMILES>
        '''
        compoundid,standard = row.strip().split("\t")

        try: self._execute(self.nsmiles.insert(), {'smiles':standard})
        except IntegrityError: pass
        standardid = self.sess.query(self.nsmiles.c.id).filter(self.nsmiles.c.smiles == standard).one()[0]
        try:
            self._execute(self.standards.insert(),
                {'compound':compoundid,
                 'standardizer':standardizerid,
                 'standard':standardid})
        except IntegrityError as ie:
            print("entry exists for compound and standardizer", compoundid, standardizerid, standardid)
        except:
            print("failed to insert ", compoundid, standardizerid, standardid)
            raise
        self.sess.commit()    

    def getsmiles(self, compid):
        return self.sess.query(self.nsmiles.c.smiles).filter(self.nsmiles.c.id == compid).one()[0]

    def getsmirks(self, reacid):
        signature = self.sess.query(self.reaction.c.signature).filter(self.reaction.c.id == reacid).one()[0]
        substrateids, productids = [ x.split('.') for x in signature.split('>>') ]
        substrates = [ self.getsmiles(x) for x in substrateids ]
        products = [ self.getsmiles(x) for x in productids ]
        return '{0}>>{1}'.format('.'.join(substrates),'.'.join(products))

    def getruleannotation(self, reacid):
        data = self._execute('''
select bt_rule
from edb_rule_reaction rr
 join edb_reaction er on er.id = rr.edb_reaction
where er.reaction = {0}
union
select container||'-'||simple
from rule_reaction x join rule on rule.id = x.rule
where reaction = {0}
'''.format(reacid))
        return [e[0] for e in data]

    def getfocus(self, reactor):
        data = self._execute('''
select rr.edb_reaction, bt_rule, r.id, r.signature,
       string_agg(cast (step_no as text),',') steps,
       string_agg(distinct container, ',') confirmed,
       count(case when reactor=1 then 1 end) envipath,
       count(case when reactor=2 then 1 end) eawag
from edb_rule_reaction rr
 join edb_reaction er on er.id = rr.edb_reaction
 join reaction r on r.id = er.reaction
 left join rule_reaction ru on ru.reaction = r.id
 left join rule on ru.rule = rule.id
group by rr.edb_reaction, bt_rule, r.id, r.signature
having count(case when reactor=2 then 1 end) = 0
order by rr.edb_reaction
''')
        return data.fetchall()

    def findequivalents(self, edb_reaction, reactor):
        data = self._execute('''
select reaction, signature, rule, container, simple, reactor
from (
select id, signature
from reaction r
 join rrole p on r.id = p.reaction and p.isproduct is true
 join rrole s on r.id = s.reaction and s.isproduct is false
where p.compound in (
   select standard 
   from rrole natural join edb_reaction natural join standards
   where id = '{0}' and isproduct is true
) and s.compound in (
   select standard 
   from rrole natural join edb_reaction natural join standards
   where id = '{0}' and isproduct is false
)) as eqmod
 join rule_reaction on eqmod.id = reaction and reactor = {1}
 join rule on rule.id = rule
'''.format(edb_reaction, reactor))
        return data.fetchall()

    def getreactions(self, reacids=None):
        """
        get reactions and signatures from reaction table
        """
        return self.sess.query(self.reaction.c.id, self.reaction.c.signature)\
            .filter(self.reaction.c.id.in_(reacids)).all()

    def getproducts(self, educts, reactorid):
        """
        get signature, products, substrate, rule  with matching reactorid and educts
        because the performance of the query decreases dramatically when 
        reactorid is an internal condition the filtering is done outside.
        """
        data = self._execute('''
select distinct signature, s.compound substrate, rule, p.compound product, reactor
from reaction join (
  select rule, reaction, reactor
  from rule_reaction
  where rule not in (select rule from rule_exclusion)
) as rr on reaction.id = rr.reaction join (
  select compound, reaction
  from rrole
  where isproduct is false
  and compound in ({0})
) as s using(reaction) join (
  select compound, reaction
  from rrole
  where isproduct is true
) as p using(reaction)
'''.format(",".join([str(e) for e in educts])))
        return [f for f in data.fetchall() if str(f['reactor']) == str(reactorid)] 

    def getstandards(self, compounds, standardizers=None):
        """
        get compound, standardizer, standard  by compounds and standardizer
        """
        if standardizers is None:
            return self.sess.query(self.standards.c.compound,
                                   self.standards.c.standardizer,
                                   self.standards.c.standard)\
                .filter(self.standards.c.compound.in_(compounds))\
                .all()
        else:
            return self.sess.query(self.standards.c.compound,
                                   self.standards.c.standardizer,
                                   self.standards.c.standard)\
                .filter(self.standards.c.compound.in_(compounds))\
                .filter(self.standards.c.standardizer.in_(standardizers))\
                .all()

    def getreversestandards(self, standards, standardizers=None):
        """
        get compound, standardizer, standard  by compounds and standardizer
        """
        if standardizers is None:
            return self.sess.query(self.standards.c.compound,
                                   self.standards.c.standardizer,
                                   self.standards.c.standard)\
                .filter(self.standards.c.standard.in_(standards))\
                .all()
        else:
            return self.sess.query(self.standards.c.compound,
                                   self.standards.c.standardizer,
                                   self.standards.c.standard)\
                .filter(self.standards.c.standard.in_(standards))\
                .filter(self.standards.standardizers.in_(standardizers))\
                .all()

    def getirrelevantcompounds(self):
        """
        compounds that are not considered relevant for a perfect match
        :return: list of compound ids
        """
        return [ x[0] for x in self.sess.query(self.irrelevant.c.id).all() ]

    def declareirrelevant(self, compound):
        """
        setter for irrelevant compounds
        :param compound:
        :return:
        """
        self._execute(self.irrelevant.insert(), {'id': compound})
        self.sess.commit()

    def excludereaction(self, edb_reaction, reason):
        """
        setter for excluded reactions
        :param edb_reaction: shortkey (reacid) of package reaction
        :param reason: catchy string, up to 8 characters
        :return:
        """
        try:
            self._execute(self.edb_reaction_exclusion.insert(), {'edb_reaction': edb_reaction, 'reason': reason})
        except IntegrityError:
            print("failed to insert, reaction {0} is already excluded".format(edb_reaction))
        self.sess.commit()

    def excluderulecontainer(self, container, reason):
        """
        setter for excluded rules
        :param container: shortkey (bt) of package container rule
        :param reason: catchy string, up to 8 characters
        :return:
        """
        for rule in self.sess.query(self.rule.c.id).filter(self.rule.c.container == container).all():
            try:
                self._execute(self.rule_exclusion.insert(), {'rule': rule[0], 'reason': reason})
            except IntegrityError:
                print("failed to insert, rule {0} from {1}, it is already excluded".format(rule[0], container))
            self.sess.commit()

    def matchreport(self, reactorid, rulepathset):
        try:
            mid = None
            for product, data in rulepathset:
                if mid is None:
                    mid = self._execute(self.match_es.insert().returning(self.match_es.c.id),
                                        {'reactor': reactorid, 'reaction': data.reaction}).fetchone().id
                mpid = self._execute(self.match_parts.insert().returning(self.match_parts.c.id),
                                     {'match': mid,
                                      'final_compound': product,
                                      'initial_compound': data.initial_compound,
                                      'initial_rule': data.initial_rule,
                                      'steps': len(data.rules),
                                      'standardized': len(data.standardizers) > 0}).fetchone().id
                i=0
                for rule in reversed(data.rules):
                    i+=1
                    self._execute(self.match_rules.insert({'match_part': mpid, 'rule': rule, 'step': i}))
                for standardizer in data.standardizers:
                    self._execute(self.match_standardizers.insert({'match_part': mpid, 'standardizer': standardizer}))
        except Exception as e:
            print("ERROR on matchreport: {0} - {1}".format(e.__class__.__name__, str(e)))
            raise
        self.sess.commit()

    def generationreport(self, reactor, reaction, generation, yieldn):
        try:
            self._execute(self.generationsize.insert({'reactor':reactor,'reaction':reaction,'generation':generation,'nproducts':yieldn}))
        except Exception as e:
            print("ERROR on generationreport: {0} - {1}".format(e.__class__.__name__, str(e)))
            raise
        self.sess.commit()

    def reactionmap(self):
        """
        :return: edb_reactions as a map with 'reacid' as key and reaction.id as value
        """
        rmap = {}
        for (reacid, reaction) in self.sess.query(self.edb_reaction.c.id, self.edb_reaction.c.reaction).all():
            rmap[reacid] = reaction
        return rmap

    def getreactorid(self, reactorname):
        return self.sess.query(self.reactor.c.id).filter(self.reactor.c.name == reactorname).one()[0]

