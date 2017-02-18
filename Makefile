#!/bin/bash

build: filtru

filtru: filter.c
	mpicc -o filtru filter.c

clean:
	rm filtru			