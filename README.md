# postgresql-signal-similarity

PostgreSQL C function for evaluating similarity between multivariate signals. Based on DTW algoritm described here http://www.cs.ucr.edu/~eamonn/UCRsuite.html .

## Synopsis

`signal_similarity(text[] query_series_names, float8[][] query_values, text[] data_series_names, float8[][] data_values, float4 warping_window)`

Return array with structure as follow:

`{similarity, query distance for series 1, 2..., location (index) of matched query on series 1, 2...}`

## Installation on linux

- install `postgres-server-dev-X.Y`


- `make` (build shared library)
- `sudo make install` (copy the new shared object (.so) file to your pgsql library directory `pg_config --pkglibdir`)
- `sudo -u postgres psql -d YOURDATABASE -a -f signal_similarity.sql `  (create function in your desired database)


## Author

Rafael Pawlos

Based on DTW algorithm described in paper:

Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen, Gustavo Batista, Brandon Westover, Qiang Zhu, Jesin Zakaria, Eamonn Keogh (2012). Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping SIGKDD 2012.