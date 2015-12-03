#!/bin/bash

id=$1

for i in $(find /data/storage/patients -name *$id* -type d)
do
	path=$(echo $i | sed 's/\/data\/storage\/patients\///')
	s3cmd sync --no-check-md5 --continue-put --recursive $i s3://averapatients/$id/$path/
done

exit 0



