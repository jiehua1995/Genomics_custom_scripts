#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
import multiprocessing
import time
import math
import platform
import argparse

# -----------------------------
# 检查依赖
# -----------------------------
required_packages = ["psutil", "tqdm"]
missing_packages = []

for pkg in required_packages:
    try:
        __import__(pkg)
    except ImportError:
        missing_packages.append(pkg)

if missing_packages:
    print("[WARNING] The following Python packages are missing: {}".format(", ".join(missing_packages)))
    print("You can install them using pip or conda:")
    for pkg in missing_packages:
        print("  pip install {}".format(pkg))
        print("  conda install {}".format(pkg))
    sys.exit(1)

# -----------------------------
# 导入依赖
# -----------------------------
import psutil
from tqdm import tqdm

# -----------------------------
# 硬件检测
# -----------------------------
def detect_hardware():
    print("\n=== Hardware Info ===")
    # CPU
    cpu_count = multiprocessing.cpu_count()
    try:
        cpu_freq = psutil.cpu_freq()
        cpu_current = cpu_freq.current
        cpu_max = cpu_freq.max
    except:
        cpu_current = cpu_max = "Unknown"
    cpu_model = platform.processor()
    print("CPU: {}".format(cpu_model))
    print("Cores/Threads: {}".format(cpu_count))
    print("Frequency: {} MHz (max: {} MHz)".format(cpu_current, cpu_max))

    # 内存
    mem = psutil.virtual_memory()
    print("Memory: {:.2f} GB".format(mem.total / 1024**3))

    # GPU (NVIDIA)
    try:
        gpus = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"]
        ).decode().strip()
        if gpus:
            print("GPU(s):")
            for line in gpus.split("\n"):
                print("  " + line)
    except:
        print("GPU: Not detected or no NVIDIA GPU")

    # 磁盘
    disk = psutil.disk_usage('/')
    print("Disk: {:.2f} GB".format(disk.total / 1024**3))
    print("====================\n")

# -----------------------------
# CPU 测试函数
# -----------------------------
def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def cpu_task(start, end):
    count = 0
    for i in range(start, end):
        if is_prime(i):
            count += 1
    return count

def cpu_task_with_progress(start, end):
    count = 0
    for i in tqdm(range(start, end), ncols=80):
        if is_prime(i):
            count += 1
    return count

# -----------------------------
# 基准测试
# -----------------------------
def run_single_thread(n):
    print("[Benchmark] Single-thread: calculating primes 0 ~ {}".format(n))
    start_time = time.time()
    count = cpu_task(0, n)
    elapsed = time.time() - start_time
    print("Found {} primes in {:.2f}s\n".format(count, elapsed))
    return elapsed

def run_multi_thread(n, threads):
    print("[Benchmark] Multi-thread ({} threads): calculating primes 0 ~ {}".format(threads, n))
    start_time = time.time()
    pool = multiprocessing.Pool(threads)
    chunk_size = n // threads
    ranges = [(i*chunk_size, (i+1)*chunk_size) for i in range(threads)]
    results = pool.starmap(cpu_task_with_progress, ranges)
    pool.close()
    pool.join()
    total_primes = sum(results)
    elapsed = time.time() - start_time
    print("Found {} primes in {:.2f}s\n".format(total_primes, elapsed))
    return elapsed

# -----------------------------
# 主程序
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Server Benchmark Script: Detect hardware and measure CPU performance (prime calculation)."
    )
    parser.add_argument(
        "-n", "--number", type=int, default=50000,
        help="Upper limit of prime calculation (default: 50000)"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=2,
        help="Number of threads for multi-thread mode (default: 2)"
    )
    parser.add_argument(
        "--allcores", action="store_true",
        help="Use all available CPU cores for multi-thread test"
    )

    args = parser.parse_args()

    detect_hardware()

    # 单线程
    run_single_thread(args.number)

    # 多线程
    thread_num = multiprocessing.cpu_count() if args.allcores else args.threads
    run_multi_thread(args.number, thread_num)

if __name__ == "__main__":
    main()
