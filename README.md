# GAWM: Genome Analyst with Mongodb

## Start mongod

```bash
# rm ~/share/mongodb/data/mongod.lock
numactl --interleave=all ~/share/mongodb/bin/mongod --config ~/share/mongodb/mongod.cnf
```

## GC wave and bed count

```bash
cd ~/Scripts/gawm

# drop db if exists
mongo S288c_gc --eval "db.dropDatabase();"

# generate
perl gen_mg.pl -d S288c_gc -n S288c --dir ~/data/alignment/example/scer/Genomes/S288c --parallel 1

# GC
perl insert_gcwave.pl -d S288c_gc --batch 1 --parallel 8

# CV
perl update_sw_cv.pl -d S288c_gc --batch 1 --parallel 8

# bed count
perl count_bed.pl -d S288c_gc --run insert -f doc/spo11_hot.bed --batch 1 --parallel 1
perl count_bed.pl -d S288c_gc --run count --batch 1 --parallel 8

# stats
perl stat_mg.pl -d S288c_gc --index
```

## genome features: ofgsw

```bash
cd ~/Scripts/gawm

mongo S288c_spo11 --eval "db.dropDatabase();"

perl gen_mg.pl -d S288c_spo11 -n S288c --dir ~/data/alignment/example/scer/Genomes/S288c --parallel 1

perl insert_bed.pl -d S288c_spo11 -tag spo11 -f doc/spo11_hot.bed --batch 1 --parallel 8

perl update_sw_cv.pl -d S288c_spo11 --batch 1 --parallel 8

perl stat_mg.pl -d S288c_spo11 --index --by tag --replace ofg="DSBs"
```

## Prof GC wave

```bash
cd ~/Scripts/gawm

# drop db if exists
mongo S288c_prof --eval "db.dropDatabase();"

# generate
perl gen_mg.pl -d S288c_prof -n S288c --dir ~/data/alignment/example/scer/Genomes/S288c/VII.fa --parallel 1

# GC
perl -d:NYTProf insert_gcwave.pl -d S288c_prof --batch 1 --parallel 4
nytprofhtml --open

# CV
perl -d:NYTProf update_sw_cv.pl -d S288c_prof --batch 1 --parallel 4
nytprofhtml --open

# stats
perl stat_mg.pl -d S288c_prof
```
