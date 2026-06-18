#### CURVE d=0 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 1.0816 | 0.0655 | 12.1897 |
| 10000 | 4.8599 | 0.5679 | 124.8448 |
| 100000 | 73.7851 | 5.6363 | 1305.2842 |
| 1000000 | 950.6003 | 56.1480 | nan |

#### CURVE d=0 — PER-CALL python loop (median ms over N calls)

| N | pygbs | scipy |
|---|---|---|
| 1000 | 0.5398 | 2.4554 |
| 10000 | 5.4481 | 24.2488 |

#### CURVE d=1 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 0.5493 | 0.0637 | nan |
| 10000 | 4.8311 | 0.5842 | nan |
| 100000 | 67.5423 | 5.8182 | nan |
| 1000000 | 962.4100 | 58.2721 | nan |

#### CURVE d=2 — BATCHED (median ms, lower=better)

| N | pygbs | scipy | geomdl |
|---|---|---|---|
| 1000 | 0.5026 | 0.0633 | nan |
| 10000 | 4.9175 | 0.5930 | nan |
| 100000 | 66.8914 | 5.9198 | nan |
| 1000000 | 964.4146 | 59.1272 | nan |
