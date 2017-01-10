# Benchmarks on different versions of MongoDB

## Test command lines

```bash
export GAWM_MONGO_DIR=$HOME/share/mongodb30/bin
export GAWM_PORT=27030
export GAWM_PARALLEL=4

cd ~/Scripts/gawm/

#----------------------------#
# gcwave
#----------------------------#
$GAWM_MONGO_DIR/mongo Atha_GC --port $GAWM_PORT --eval "db.dropDatabase();"

# 257: db.getCollection('align').count({})
perl gen_mg.pl -d Atha_GC -n Atha --port $GAWM_PORT \
    --dir ~/data/alignment/Ensembl/Atha/ \
    --length 500000 --parallel $GAWM_PARALLEL

perl insert_gcwave.pl -d Atha_GC --port $GAWM_PORT --batch 10 --parallel $GAWM_PARALLEL

perl update_sw_cv.pl -d Atha_GC --port $GAWM_PORT --batch 10 --parallel $GAWM_PARALLEL

perl stat_mg.pl -d Atha_GC --port $GAWM_PORT --chart

#----------------------------#
# T-DNA ofgsw
#----------------------------#
$GAWM_MONGO_DIR/mongo Atha_TDNA_SW --port $GAWM_PORT --eval "db.dropDatabase();"

perl gen_mg.pl -d Atha_TDNA_SW -n Atha --port $GAWM_PORT \
    --dir ~/data/alignment/Ensembl/Atha/ \
    --length 500000 --parallel $GAWM_PARALLEL

perl insert_position.pl -d Atha_TDNA_SW --port $GAWM_PORT --batch 10 --parallel $GAWM_PARALLEL \
    --style center \
    --tag tdna --type CSHL -f ~/data/salk/Atha/T-DNA.CSHL.pos.txt \
    --tag tdna --type FLAG -f ~/data/salk/Atha/T-DNA.FLAG.pos.txt \
    --tag tdna --type MX   -f ~/data/salk/Atha/T-DNA.MX.pos.txt   \
    --tag tdna --type RATM -f ~/data/salk/Atha/T-DNA.RATM.pos.txt

perl stat_mg.pl -d Atha_TDNA_SW --port $GAWM_PORT --by type --chart --replace ofg="insert sites"

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

* 4 threads: i7-6700k, 32G, SSD, macOS 10.11

    |        |          |  2.6.11 |   3.0.7 |   3.4.1 |
    |:------:|:--------:|--------:|--------:|--------:|
    | gcwave |   gen    |     9'' |     9'' |     9'' |
    |        |  gcwave  | 10'52'' | 10'53'' | 10'48'' |
    |        |  sw_cv   | 41'19'' | 40'55'' | 40'56'' |
    |        |   stat   |         |    12'' |    10'' |
    | ofgsw  | position |  5'16'' |  4'58'' |  4'45'' |
    |        |   stat   |    21'' |    27'' |    22'' |
    |        |          |         |         |         |
    |        |          |         |         |         |
    |        |          |         |         |         |

