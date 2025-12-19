#!/bin/bash

# ===============================
# Настройки эксперимента
# ===============================

THREADS=(1 2 4 8 16)
TOTAL_OPS=100000
INSERTS_IN_MAIN=1000
RUNS=5

# Сценарии нагрузки
# формат: name search insert
SCENARIOS=(
  "read_heavy 0.9 0.05"
  "balanced   0.5 0.25"
  "write_heavy 0.2 0.4"
)

# Выходной CSV
OUTFILE="results.csv"

# ===============================
# Компиляция
# ===============================

echo "Compiling..."
gcc -O2 -std=c11 -pthread pth_ll_rwl.c my_rand.c -o pthread_rwlock
gcc -O2 -std=c11 -pthread pth_ll_my_rwl.c my_rand.c -o my_rwlock

if [ $? -ne 0 ]; then
    echo "Compilation failed"
    exit 1
fi

# ===============================
# Заголовок CSV
# ===============================

echo "impl,threads,scenario,run,time,member,insert,delete" > $OUTFILE

# ===============================
# Запуск тестов
# ===============================

for impl in pthread my; do
  if [ "$impl" == "pthread" ]; then
    BIN=./pthread_rwlock
  else
    BIN=./my_rwlock
  fi

  for t in "${THREADS[@]}"; do
    for sc in "${SCENARIOS[@]}"; do
      set -- $sc
      SC_NAME=$1
      SEARCH=$2
      INSERT=$3

      for ((r=1; r<=RUNS; r++)); do
        echo "Running $impl | threads=$t | $SC_NAME | run=$r"

        OUTPUT=$(printf "%d\n%d\n%f\n%f\n" \
          $INSERTS_IN_MAIN $TOTAL_OPS $SEARCH $INSERT | \
          $BIN $t)

        TIME=$(echo "$OUTPUT" | grep "Elapsed time" | awk '{print $4}')
        MEMBER=$(echo "$OUTPUT" | grep "member ops" | awk '{print $4}')
        INS=$(echo "$OUTPUT" | grep "insert ops" | awk '{print $4}')
        DEL=$(echo "$OUTPUT" | grep "delete ops" | awk '{print $4}')

        echo "$impl,$t,$SC_NAME,$r,$TIME,$MEMBER,$INS,$DEL" >> $OUTFILE
      done
    done
  done
done

echo "Done. Results saved to $OUTFILE"
