#### CURVE d=0 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 0.1397 | 0.0599 | 12.3131 |
| 10000 | 1.3984 | 0.5863 | 120.5931 |
| 100000 | 9.0879 | 5.5896 | 1264.0852 |
| 1000000 | 92.0227 | 56.0400 | nan |

#### CURVE d=0 — PER-CALL python loop (median ms over N calls)

| N | pygbs | scipy |
|---|---|---|
| 1000 | 0.5459 | 2.4143 |
| 10000 | 5.3642 | 24.1735 |

#### CURVE d=1 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 0.1348 | 0.0613 | nan |
| 10000 | 1.0366 | 0.6005 | nan |
| 100000 | 9.1948 | 5.7927 | nan |
| 1000000 | 101.2452 | 58.2816 | nan |

#### CURVE d=2 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 0.1377 | 0.0616 | nan |
| 10000 | 1.0641 | 0.6071 | nan |
| 100000 | 9.4953 | 5.8908 | nan |
| 1000000 | 100.4459 | 58.9033 | nan |
