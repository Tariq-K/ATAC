#####################################################
####              Helper functions               ####
#####################################################

import cgatcore.pipeline as P
import glob
import pandas as pd
import sqlite3
import os

# Pipeline configuration
P.get_parameters(
		 ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
		  "../pipeline.yml",
		  "pipeline.yml"],
		 )

PARAMS = P.PARAMS

db = PARAMS['database']['url'].split('./')[1]

def connect():
    '''connect to database.
    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(db)

    if not os.path.exists(PARAMS["annotations_database"]):
        raise ValueError(
                     "can't find database '%s'" %
                     PARAMS["annotations_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
    (PARAMS["annotations_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh



def fetch_DataFrame(query,
                    dbhandle=db):
    '''Fetch query results and returns them as a pandas dataframe'''

    dbhandle = sqlite3.connect(dbhandle)

    cc = dbhandle.cursor()
    sqlresult = cc.execute(query).fetchall()
    cc.close()

    # see http://pandas.pydata.org/pandas-docs/dev/generated/
    # pandas.DataFrame.from_records.html#pandas.DataFrame.from_records
    # this method is design to handle sql_records with proper type
    # conversion

    field_names = [d[0] for d in cc.description]
    pandas_DataFrame = pd.DataFrame.from_records(
        sqlresult,
        columns=field_names)
    return pandas_DataFrame


# ---------------------------------------------------
