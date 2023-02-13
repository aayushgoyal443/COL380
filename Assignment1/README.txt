I have used different approaches to optimize the result. I will quickly go through them.

Approach1
1:Idea
-- The very basic approach that I used was to construct the entire matrix and then mutiply it with itself. Now this approach is correct but it is very space ineefiecient as well as time inefficient. The space complexity is O(n^2) and the time complexity is O(n^3). 

1:Why Implemented
-- I implemented this approach to get a basic understanding of the problem and to get a feel of the problem. Also using this tried how to use pragma omp loops and tasts in the C++ code. So this was a way for me to get familiar with the OpenMP library.

1:Result 
-- Result was it was able to multiply the matrix in input1 in around 20ms on my system. But for the input2 I have no stats since process was getting killed in between. Hence I had to abandon this approach.

1:Drawback
-- The main drawback of this is that it requires us to construct the entire matrix and then multiply them normally. This defeats the entire purpose of giving enput as a sparse matrix and also asking for the output as a sparse matrix. Moreover this solution is not scalable as the space complexity is O(n^2) and the time complexity is O(n^3). For the input2 where the value of n is 10**5, the space complexity is 10^10 and the time complexity is 10^15. This is not feasible for any practical application and I am not getting the output either. Hence I changed my approach 

Approach 2 and the Final Approach
2: Idea
-- The idea here is to multiply them like sparse matrices without converting them into the original matrix form. So for this I made a class called Block and another class called Sparse Matrix. IN the sparse matrix I store i,j, and reference to the block. The block class basically stores the index where the block starts in the matrix and then the block. For all the indices where i!=j, I have also stored it's transpose matrix. I have stored it already because that makes the work easy for us when we are doing the main part of multiplication. Doing these transpose and checks will take a lot of time if done multiple times during the multiplication, and hence it is better to store them once and for all.  For the multiplication of blocks I have kept the same approach as the previous approach, because I am assuming the size of the blocks will be small. 
Now coming to the main part about the multiplication of sparse matrix with itself, I started travesing through the list of blocks in one loop. Now for each block I traverse again through the list of blocks and check if it is possible to multiply them. Now if it possible to multiply them, does it belong to the upper traingular region only? If yes then I have multipled the blocks and store the result in the output desired matrix.

For Parallelizing:
Now since we are travesing through the list of blocks, this can be easily divided into different task. I can individual task for each pair of block*block. Now these tasks will get scheduled on the different threads that are available in the system. In the end they all wil combine and give us the final result. In between that we have to take care about the cases of shared memory. For that I have used critical.

Critical section:
I have also used "#pragma omp critical (name)" to make sure that the output matrix is not accessed by multiple threads at the same time. Critical section is basically around the point where multiple threads can be trying to write to the same block of the output matrix. So I have used this to make sure that only one thread is writing to a particular block of the output matrix at a time. I kept name as "i_j" where i and j are the indices of the block. This way I can make sure that different locks are there for every block. Hence different blocks of the same resultant matrix can be written at the same time.


2: Why Implemented
-- I implemented this approach because it ewas easy to implement as well as scalable. Since in this approach we are not converting the sparse matrix into the original matrix form, we are not wasting any space. Also we are not doing any extra work. We are just multiplying the blocks and storing the result in the output matrix which is also a sparse matrix. Space complexity will be O(m^2*k^2) and time complexity will be O(m^3*k^3). Since m is supposed to be small this scalable solution.


2: Result
My code is working for both the input1 and the input2 now. And it is giving output in reasonable time. For input1 it is taking around 1ms and for input2 it is taking around on the css machines. Before doing parallelization my code was taking around 200s on the input2 but after parallelization it is taking around 17s only which is a very good improvement. So the parallelization has helped a lot in making it faster. Space complexity is also as efficient as possible, I could not think of any ways to recduce the space complexity further.


2: Drawback
A drawback I can still think is that I have not parallelized the part where I am reading input and storing their transpose matrix. This might or might not reduce the runtime, since reading input and storing them in class is already taking negligible time compared to the multiplication part, but due to overhead of thread management, it is sometimes taking more time than originally. Hence I didn't parallelize this part.
Another drawback could be I am not using the O(n^(log_sub{2}_sup{7})) algorithm for multiplication of blocks. But I have not done this because the value of m is very small and then this adds unnecessary reccursive calls and we know that there is again some over head of doing a lot of recursive function calls. Hence I went for the normal O(n**3) approach for multiplication of blocks.


Final scalablitity Analysis:


