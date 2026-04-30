import csv
import gzip
import sys

import pyarrow.parquet as pq

src = sys.argv[1]
dst = sys.argv[2]
mode = sys.argv[3]
pf = pq.ParquetFile(src)
schema = pf.schema_arrow

if mode == "parquet":
    writer = pq.ParquetWriter(dst, schema, compression="snappy", use_dictionary=True)
    try:
        for batch in pf.iter_batches(batch_size=65536):
            writer.write_batch(batch)
    finally:
        writer.close()
elif mode == "csv_gz":
    with gzip.open(dst, "wt", newline="") as fh:
        writer = csv.writer(fh)
        names = schema.names
        writer.writerow(names)
        for batch in pf.iter_batches(batch_size=65536):
            cols = [batch.column(i).to_pylist() for i in range(batch.num_columns)]
            if not cols:
                continue
            for row in zip(*cols):
                writer.writerow(row)
else:
    raise SystemExit(f"unsupported mode: {mode}")
