#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <unistd.h>
#include <thread>
#include <omp.h>
using namespace std;
using namespace std::chrono;

#define MATSIZE 100


//this calculates the total of the matrix to confirm that we arent doing any improper maths and that the rounding error is close enough (with doubles there is no rounding error)
void calculate_total(double *array)
{
    double total = 0;
    for(int i = MATSIZE-1; i>=0; i--) //backwards loop
    {
            for (int j = MATSIZE -1 ; j>=0; j--) //backwards loop
            {
                total += array[i*MATSIZE + j];
            }
    }

    cout << "Total is: " << total << endl;
}


//this prints the matrix with colored backgrounds, doesn't work too well on large arrays because it takes like three lines per line of the matrix
void print_matrix(double *array)
{
    cout << endl;
    for (int i = 0; i<MATSIZE; i++){
        cout << "[";
        for (int j =0; j<MATSIZE; j++){
            double val = array[i * MATSIZE + j] ; 
            if (val >= 100){
                cout << "\033[1;41m " << val << "\033[m\t";
            }
            else if (val > 50){
                cout << "\033[1;42m " << val << "\033[m\t";
            }
            else if (val > 30){
                cout << "\033[1;43m " << val << "\033[m\t";

            }
            else if (val > 20){
                cout << "\033[1;44m " << val << "\033[m\t";
            }
            else if (val > 12){
                cout << "\033[1;45m " << val << "\033[m\t";
            }
            else if (val > 8){
                cout << "\033[1;46m " << val << "\033[m\t";
            }
            else if (val > 1){
                cout << "\033[1;47;30m " << val << "\033[m\t";
            }
            else {
                cout << "\033[0m " << val << "\033[m\t";
            }
        }
        cout << "]\n" ;
    }
    cout << endl;
}

// this prints the matrix as different colored spaces with different colors representing different value ranges
void print_matrix_colors(double *array)
{
    cout << endl;
    for (int i = 0; i<MATSIZE; i++){
        cout << "[";
        for (int j =0; j<MATSIZE; j++){
            double val = array[i * MATSIZE + j] ; 
            if (val > 30){
                cout << "\033[1;41m \033[m";
            }
            else if (val > 15){
                cout << "\033[1;42m \033[m";
            }
            else if (val > 14){
                cout << "\033[1;43m \033[m";

            }
            else if (val > 13){
                cout << "\033[1;44m \033[m";
            }
            else if (val > 12){
                cout << "\033[1;45m \033[m";
            }
            else if (val > 8){
                cout << "\033[1;46m \033[m";
            }
            else if (val > 1){
                cout << "\033[1;47;30m \033[m";
            }
            else {
                cout << "\033[0m \033[m";
            }
        }
        cout << "]\n" ;
    }
    cout << "===================================" << endl;
}


// this is the main function that we call on every cell to find it and its neighbors new values
void compute_cell(double *array, int x, int y)
{
    
    double avg = 0;
    int count = 0;  
    // doing the x and y offsets this way instead of doing multiple for loops potentially saves a bit of time from avoiding the branches
    int xoffset[8] = {1,1,1,0,0,-1,-1,-1};
    int yoffset[8] = {1,0,-1,1,-1,1,0,-1};

    // look to see if we are on the border, if so we need to be careful, otherwise we know that we dont have to range check so we can save some clock cycles and take the cheaper route below
    if (x == 0 || y == 0 || x == MATSIZE-1 || y == MATSIZE-1)
    {
        // this is when we are near the border to compute the average
        for (int i = 0; i < 8; i++)
        {
            // using the offset arrays we can make this a single for loop and can skip having to check if we are in the middle of the neighbors (the cell itself) saving us another branch
            int currx = x + xoffset[i];
            int curry = y + yoffset[i];

            // some booleans that make the next if statement a whole lot cleaner
            bool Xabove = currx >= 0;
            bool XbelowMATSIZE = currx < MATSIZE;
            bool Yabove = curry >= 0;
            bool YbelowMATSIZE = curry < MATSIZE;

            //if we are within the bounds of the array then find the value at that array and add it to the average
            if (Xabove && XbelowMATSIZE && Yabove && YbelowMATSIZE)
            {
                avg += array[(currx)*MATSIZE + (curry)];
                count +=1;
            }
        }
        
        // this is now where we distribute the values
        // calculate what we should be sending to neighbors (assuming they are all above 0 after moving the units)
        double difference = (array[x*MATSIZE + y] - (avg/count)) * 0.05; 

        
  
        //reset the count to 0 so we can track any possible negatives
        count = 0;

        // go through all the neighbors and see if we can move units
        for (int i = 0; i < 8; i++)
        {
            // using the offset arrays we can make this a single for loop and can skip having to check if we are in the middle of the neighbors (the cell itself) saving us a branch
            int currx = x + xoffset[i];
            int curry = y + yoffset[i];

            //some bools that make the next if cleaner
            bool Xabove = currx >= 0;
            bool XbelowMATSIZE = currx < MATSIZE;
            bool Yabove = curry >= 0;
            bool YbelowMATSIZE = curry < MATSIZE;

            //if we are within the bounds of the array then move the difference
            if (Xabove && XbelowMATSIZE && Yabove && YbelowMATSIZE)
            {
                // if distributing the units wont cause the neighbor to go into the negatives then do so
                if (array[(currx)*MATSIZE + curry] + difference > 0)
                {
                    array[(currx)*MATSIZE + (curry)] += difference;
                    count +=1;
                }
            }
        }
        // calculate the new value for the current cell based on how many units we were able to move
        array[x*MATSIZE + y] = array[x*MATSIZE + y] -  (difference * count);
    }
    else 
    // when we arent near the border we dont need to check each cell to see if we are within bounds as we already know that so we can save quite a bit of bound checking
    {
        // go through all the neighbors
        for (int i = 0; i < 8; i++)
        {
            // using the offset arrays we can make this a single for loop and can skip having to check if we are in the middle of the neighbors (the cell itself) saving us a branch
            int currx = x + xoffset[i];
            int curry = y + yoffset[i];

            // there is no need to check to see if we are within bounds here so this will save quite a few instructions
            avg += array[(currx)*MATSIZE + (curry)];
            count +=1;
        }
        
        // this is now where we distribute the values
        double difference = (array[x*MATSIZE + y] - (avg/count)) * 0.05; 

        //reset the count to 0 so we can count the number of non 0 totals
        count = 0;
        for (int i = 0; i < 8; i++)
        {
            // using the offset arrays we can make this a single for loop and can skip having to check if we are in the middle of the neighbors (the cell itself) saving us a branch
            int currx = x + xoffset[i];
            int curry = y + yoffset[i];

            if (array[(currx)*MATSIZE + curry] + difference > 0)
            {
                array[(currx)*MATSIZE + (curry)] += difference;
                count +=1;
            }
        }
        // take the units away from this cell (or add them if it pulled units)
        array[x*MATSIZE + y] = array[x*MATSIZE + y] -  (difference * count);
    }
}

// this is the single threaded implementation
double do_st()
{
    //get the current time
    auto starttime = chrono::steady_clock::now();

    // we create an array in a contiguous block of memory as a 1d array but when accessing it instead of [i][j] we do [i*rowsize + j], this makes it much faster as then all of the array is in one contiguous block of memory
    double *array;
    array = new double[MATSIZE * MATSIZE];

    //dump 100000 into the middle of the leftmost column of the array
    int midpoint = static_cast<int>(MATSIZE/2); 
    array[midpoint*MATSIZE] = 100000;

    // set the current max to some arbitrary value above 12
    double curr_max = 20;
    double curr_val = 0;
    int num_iterations = 0;

    //this loop will run until the current max is below 12.
    // current max only gets computed every here and there to save cycles when we are going to be doing like 20million loops it doesnt really matter
    do
    {  
        num_iterations++;


        // this goes through all of the cells in the array, rather top down, left to right or the reverse of either of those (some ways are better but it doesnt super matter)
        //for (int i = 0; i < MATSIZE; i++) //forwards loop
        for(int i = MATSIZE-1; i>=0; i--) //backwards loop
        {
            for (int j = MATSIZE -1 ; j>=0; j--) //backwards loop
            //for (int j = 0; j < MATSIZE; j++) //forwards loop
            {
                // compute the new value of the cell and its neighbors
                compute_cell(array, i, j);
            }
        }

        // every 50000 iterations, compute the new max value of the cells. We cannot do this as we are computing the cells because they all change each other so that isnt a proper maximum at the end of computing all of the cells
        // we only do this every 50k iterations to save going through the matrix twice each run through
        if (num_iterations%50000==0){
            //set the max to the first element in the array because the last will probably be higher than anything else in this matrix
            curr_max = array[0];
        
            // this will go through and look for the maximum value in the array
            for (int i = 0; i < MATSIZE; i++) //forwards loop
            {
                for (int j = 0; j < MATSIZE; j++) //forwards loop
                {
                    // get the value of the current cell and compare it to the max that we have
                    curr_val = array[i*MATSIZE + j];

                    if (curr_val > curr_max){
                        curr_max = curr_val;
                    }
                }
            }


            // print out a little status update as well as what the matrix currently looks like
            print_matrix_colors(array);
            //print_matrix(array); //this prints the actual values but is messy with big arrays
            cout << "so far: " << num_iterations << " iterations" << "\n";
            cout << "the max value in the array is: " << curr_max << endl;

            //std::this_thread::sleep_for(std::chrono::milliseconds(2000)); // this lets us sleep during debugging
        }
    } while (curr_max >= 12.0);
    cout << endl << "it took a total of " << num_iterations << " iterations to get all of the values below 12" << endl;
    
    auto endtime = chrono::steady_clock::now();

    return chrono::duration_cast<chrono::milliseconds>(endtime-starttime).count();
}

double do_mt()
{
    //get the current time
    auto starttime = chrono::steady_clock::now();

    // we create an array in a contiguous block of memory as a 1d array but when accessing it instead of [i][j] we do [i*rowsize + j], this makes it much faster as then all of the array is in one contiguous block of memory
    double *array;
    array = new double[MATSIZE * MATSIZE];

    //dump 100000 into the middle of the leftmost column of the array
    int midpoint = static_cast<int>(MATSIZE/2); 
    array[midpoint*MATSIZE] = 100000;

    // set the current max to some arbitrary value above 12
    double curr_max = 20;
    double curr_val = 0;
    int num_iterations = 0;
    int flag = 1;

    //this loop will run until the current max is below 12.
    // current max only gets computed every here and there to save cycles when we are going to be doing like 20million loops it doesnt really matter
    do
    {  
        num_iterations++;


        // this goes through all of the cells in the array, rather top down, left to right or the reverse of either of those (some ways are better but it doesnt super matter)
        //for (int i = 0; i < MATSIZE; i++) //forwards loop
        #pragma omp parallel for
        for(int i = MATSIZE-1; i>=0; i--) //backwards loop
        {
            for (int j = MATSIZE -1 ; j>=0; j--) //backwards loop
            //for (int j = 0; j < MATSIZE; j++) //forwards loop
            {
                // compute the new value of the cell and its neighbors
                compute_cell(array, i, j);
            }
            // we need to do this to enforce the threads to meet up at the end of each row so that we dont have one thread working on row 10 and the next working on row 11 causing some race conditions
            // there is maybe a better way but the bandaid approach gets the job done
            #pragma omp critical
            flag = 1;
        }

        // every 50000 iterations, compute the new max value of the cells. We cannot do this as we are computing the cells because they all change each other so that isnt a proper maximum at the end of computing all of the cells
        // we only do this every 50k iterations to save going through the matrix twice each run through
        if (num_iterations%50000==0){
            //set the max to the first element in the array because the last will probably be higher than anything else in this matrix
            curr_max = array[0];
        
            // this will go through and look for the maximum value in the array
            // we dont need to worry about race conditions here because we arent changing anything
            #pragma omp parallel for reduction(max:curr_max)
            for (int i = 0; i < MATSIZE; i++) //forwards loop
            {
                for (int j = 0; j < MATSIZE; j++) //forwards loop
                {
                    // get the value of the current cell and compare it to the max that we have
                    curr_val = array[i*MATSIZE + j];

                    if (curr_val > curr_max){
                        curr_max = curr_val;
                    }
                }
            }


            // print out a little status update as well as what the matrix currently looks like
            print_matrix_colors(array);
            //print_matrix(array); //this prints the actual values but is messy with big arrays
            cout << "so far: " << num_iterations << " iterations" << "\n";
            cout << "the max value in the array is: " << curr_max << endl;

            //std::this_thread::sleep_for(std::chrono::milliseconds(2000)); // this lets us sleep during debugging
        }
    } while (curr_max >= 12.0);
    cout << endl << "it took a total of " << num_iterations << " iterations to get all of the values below 12" << endl;
    
    auto endtime = chrono::steady_clock::now();

    return chrono::duration_cast<chrono::milliseconds>(endtime-starttime).count();
}

int main()
{   
    // this makes sure that we always print out exactly 2 decimal places
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    
    //do this single threaded
    //double time = do_st();

    //do this multi threaded
    double time = do_mt();



    cout << "That took: " << time << "ms" << endl;
    return 0;
    
}
