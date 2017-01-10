# Benchmarks on different versions of MongoDB

## gcwave

```bash
export GAWM_MONGO_DIR=$HOME/share/mongodb/bin
export GAWM_PORT=27017
export GAWM_PARALLEL=4

cd ~/Scripts/gawm/

$GAWM_MONGO_DIR/mongo Atha --port $GAWM_PORT --eval "db.dropDatabase();"

# 257: db.getCollection('align').count({})
perl gen_mg.pl -d Atha -n Atha --port $GAWM_PORT \
    --dir ~/data/alignment/Ensembl/Atha/ \
    --length 500000 --parallel $GAWM_PARALLEL

perl insert_gcwave.pl -d Atha --port $GAWM_PORT --batch 10 --parallel $GAWM_PARALLEL

perl update_sw_cv.pl -d Atha --port $GAWM_PORT --batch 10 --parallel $GAWM_PARALLEL

perl stat_mg.pl -d Atha --port $GAWM_PORT

unset GAWM_MONGO_DIR
unset GAWM_PORT
unset GAWM_PARALLEL
```

## Start MongoDB services

```bash
# mongodb26
rm ~/share/mongodb26/data/mongod.lock
~/share/mongodb26/bin/mongod --config ~/share/mongodb26/mongod.cnf

# mongodb30
rm ~/share/mongodb30/data/mongod.lock
~/share/mongodb30/bin/mongod --config ~/share/mongodb30/mongod.cnf

# mongodb34
rm ~/share/mongodb/data/mongod.lock
~/share/mongodb/bin/mongod --config ~/share/mongodb/mongod.cnf
```

## Results

i7-6700k, 32G, SSD, macOS 10.11

|        |        | 2.6.11                    | 3.0.7                     | 3.4.1                     |
|:------:|:------:|:--------------------------|:--------------------------|:--------------------------|
| gcwave |  gen   | 9 seconds                 | 9 seconds                 | 9 seconds                 |
|        | gcwave | 10 minutes and 52 seconds | 10 minutes and 53 seconds | 10 minutes and 48 seconds |
|        | sw_cv  | 41 minutes and 19 seconds | 40 minutes and 55 seconds | 40 minutes and 56 seconds |
|        |  stat  |                           | 12 seconds                | 10 seconds                |
|        |        |                           |                           |                           |
|        |        |                           |                           |                           |
