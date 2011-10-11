#!/bin/bash 

find . -name \*.pyc -exec rm {} \;
find . -name \*~ -exec rm -i {} \;
find . -name \*.class -exec rm -i {} \;

chmod -R og-w  .
