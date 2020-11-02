import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq
import numpy as np
import os

outdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer/parquet_test'

df = pd.DataFrame({"Day":{"0":1,"1":2,"2":3,"3":4,"4":5,"5":6,"6":7,"7":8},"Year":{"0":2018,"1":2018,"2":2018,"3":2018,"4":2018,"5":2018,"6":2018,"7":2018},"Month":{"0":1,"1":1,"2":1,"3":1,"4":1,"5":1,"6":1,"7":1},"randNumCol":{"0":2,"1":5,"2":4,"3":3,"4":3,"5":5,"6":4,"7":3},"uuid":{"0":"456578af-8953-4cf7-ac27-70309353b72c","1":"df6a30da-619e-4594-a051-4fdb3572eb49","2":"7cfe724a-a827-47b1-a691-c741f4f1101d","3":"f1796ed1-f7ce-4b49-ba64-6aacdca02c0a","4":"827e4aae-1214-4c0f-ac7f-9439e8a577af","5":"08dc3c2b-b75c-4ac6-8a38-0a44007fdeaf","6":"54f4e7bb-6fd8-4913-a2c3-69ebc13dc9a2","7":"eda1dbfe-ad08-4067-b064-bcc689fa0225"},"Date":{"0":1514764800000,"1":1514764800000,"2":1514764800000,"3":1514764800000,"4":1514764800000,"5":1514764800000,"6":1514764800000,"7":1514764800000}})
table = pa.Table.from_pandas(df)
pq.write_to_dataset(table,root_path=os.path.join(outdir, 'output.parquet'),partition_cols=[ 'Day', 'Year', 'Month', 'randNumCol', 'uuid', 'Date'])

col = dict({"0":1514764800000,"1":1514851200000,"2":1514937600000,"3":1515024000000,"4":1515110400000,"5":1515196800000,"6":1515283200000,"7":1515369600000})
df = pd.DataFrame({ "NEWCOLUMN": col, "Day":{"0":1,"1":2,"2":3,"3":4,"4":5,"5":6,"6":7,"7":8},"Year":{"0":2018,"1":2018,"2":2018,"3":2018,"4":2018,"5":2018,"6":2018,"7":2018},"Month":{"0":1,"1":1,"2":1,"3":1,"4":1,"5":1,"6":1,"7":1},"randNumCol":{"0":2,"1":5,"2":4,"3":3,"4":3,"5":5,"6":4,"7":3},"uuid":{"0":"456578af-8953-4cf7-ac27-70309353b72c","1":"df6a30da-619e-4594-a051-4fdb3572eb49","2":"7cfe724a-a827-47b1-a691-c741f4f1101d","3":"f1796ed1-f7ce-4b49-ba64-6aacdca02c0a","4":"827e4aae-1214-4c0f-ac7f-9439e8a577af","5":"08dc3c2b-b75c-4ac6-8a38-0a44007fdeaf","6":"54f4e7bb-6fd8-4913-a2c3-69ebc13dc9a2","7":"eda1dbfe-ad08-4067-b064-bcc689fa0225"},"Date":{"0":1514764800000,"1":1514764800000,"2":1514764800000,"3":1514764800000,"4":1514764800000,"5":1514764800000,"6":1514764800000,"7":1514764800000}})
table = pa.Table.from_pandas(df)


pq.write_to_dataset(table,root_path=os.path.join(outdir, 'output.parquet'),partition_cols=[ 'Day', 'Year', 'Month', 'randNumCol', 'uuid', 'Date'])
