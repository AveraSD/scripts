#!/bin/bash

module load controlfreec
module load samtools

config=/data/storage/14MS10038/14MS10038_freec.txt

freec -conf $config

exit 0
