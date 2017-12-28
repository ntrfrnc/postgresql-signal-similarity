DROP FUNCTION signal_similarity(text[], float8[][], text[], float8[][], float4);

CREATE OR REPLACE function signal_similarity(text[], float8[][], text[], float8[][], float4)
RETURNS float8[]
AS 'signal_similarity', 'signal_similarity'
LANGUAGE 'c' IMMUTABLE STRICT PARALLEL SAFE
COST 500;
