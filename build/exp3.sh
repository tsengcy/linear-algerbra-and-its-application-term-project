for N in 2 10 100 1000 5000 10000; do
    printf "loop %s: " $N
    ./data $N 100
    ./LU $N 100
    ./gmres $N 100
    ./LU $N 1
    ./gmres $N 1
    ./LU $N 500
    ./gmres $N 500
    rm -f A.txt b.txt
done