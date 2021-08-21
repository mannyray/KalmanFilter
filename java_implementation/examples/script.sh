#!/bin/bash

javac test2.java ../*.java
mv ../*.class .
java test2
rm *.class
python python_verify.py
