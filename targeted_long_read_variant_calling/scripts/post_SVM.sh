#!/bin/bash

file=$1

awk '{print $1,$2}' "${file}"
