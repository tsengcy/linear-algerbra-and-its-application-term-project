for N in 2 4 5; do
    counter=0
    while [ $counter -le 100000 ]; do
        counter=`expr $counter + 1`
        ./data $N $counter
        # ./LU $N $counter
        ./Jacobi $N $counter
        # # printf "GS loop %-3s: " $counter 
        ./sor $N $counter 0.5 sor_0.5.csv
        ./sor $N $counter 1.5 sor_1.5.csv
        ./sor $N $counter 0.1 sor_0.1.csv
        ./GS $N $counter
        # ./gmres $N $counter
        rm -f A.txt b.txt
    done

done