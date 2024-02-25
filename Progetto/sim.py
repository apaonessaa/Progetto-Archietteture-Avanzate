import subprocess
import os
import time

def execute(cnt, file_name, cmd):
    with open(file_name, "w") as file:
        while cnt > 0:
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
            output = p.stdout.strip()
            file.write(f"{output}\n")
            cnt-=1

def build32(cnt, version, path, dataset, labels, k, mode):
    print(f"Start {version} ...")
    os.system("nasm -f elf32 sseutils32.nasm && nasm -f elf32 cfs32.nasm")
    os.system("gcc -m32 -msse -O0 -no-pie sseutils32.o cfs32.o "+version+".c -o "+version+" -lm")
    execute(cnt, path+version+".txt", "./"+version+" -ds "+dataset+" -labels "+labels+" -k "+k+" -"+mode+"")
    print(f"End {version} ...")

def build64(cnt, version, path, dataset, labels, k, mode):
    print(f"Start {version} ...")
    os.system("nasm -f elf64 sseutils64.nasm && nasm -f elf64 cfs64.nasm")
    os.system("gcc -m64 -msse -mavx -O0 -no-pie sseutils64.o cfs64.o "+version+".c -o "+version+" -lm")
    execute(cnt, path+version+".txt", "./"+version+" -ds "+dataset+" -labels "+labels+" -k "+k+" -"+mode+"")
    print(f"End {version} ...")

def build32omp(cnt, version, path, dataset, labels, k, mode):
    print(f"Start {version} ...")
    os.system("nasm -f elf32 sseutils32.nasm && nasm -f elf32 cfs32.nasm")
    os.system("gcc -m32 -msse -O0 -no-pie -fopenmp sseutils32.o cfs32.o "+version+".c -o "+version+" -lm")
    execute(cnt, path+version+".txt", "./"+version+" -ds "+dataset+" -labels "+labels+" -k "+k+" -"+mode+"")
    print(f"End {version} ...")

def build64omp(cnt, version, path, dataset, labels, k, mode):
    print(f"Start {version} ...")
    os.system("nasm -f elf64 sseutils64.nasm && nasm -f elf64 cfs64.nasm")
    os.system("gcc -m64 -msse -mavx -O0 -no-pie -fopenmp sseutils64.o cfs64.o "+version+".c -o "+version+" -lm")
    execute(cnt, path+version+".txt", "./"+version+" -ds "+dataset+" -labels "+labels+" -k "+k+" -"+mode+"")
    print(f"End {version} ...")

# test 5000x50
n=10
ds64="dataset/test_5000_50_64.ds"
labels64="dataset/test_5000_50_64.labels"
ds32="dataset/test_5000_50_32.ds"
labels32="dataset/test_5000_50_32.labels"
dir="./results/"
k="5"
mode="d"

build32(n, "cfs32c", dir, ds32, labels32, k, mode)

build64(n, "cfs64c", dir, ds64, labels64, k, mode)

# build32omp(n, "cfs32c_omp", dir, ds32, labels32, k, mode)

# build64omp(n, "cfs64c_omp", dir, ds64, labels64, k, mode)

# NB errore nella misurazione dei test con omp, invece che del clock() si deve usare la system per la misurazione

