### Helper functions for pipeline_motifenrichment.py
import sys
import os
import sqlite3
import cgatcore.experiment as E
from cgatcore import pipeline as P
import glob
import pandas as pd

# Pipeline configuration
P.get_parameters(
		 ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
		  "../pipeline.yml",
		  "pipeline.yml"],
		 )

PARAMS = P.PARAMS

db = PARAMS['database']['url'].split('./')[1]

def fetch(query, dbhandle=None, attach=False):
    '''Fetch all query results and return'''

    cc = dbhandle.cursor()

    if attach:
        db_execute(cc, attach)

    sqlresult = cc.execute(query).fetchall()
    cc.close()
    return sqlresult


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
