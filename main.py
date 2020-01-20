def main(Data A, B):
    # UNTITLED2 Summary of this function goes here
    # Detailed explanation goes here
    import numpy as np
    [r, c] = Data.shape;
    Score = numpy.zeros((c,c),dtype=numpy.float64);
    for i in range(1,c):
        for j in range((i + 1),c):        # Prepare the Table to find indiscarnibility%
            Table[:, 1] = Data[, i]
	    Table[:, 2] = Data[, j]
            
    print 'Hold on.. I am busy in calculating Score for i='i,'j='j
            #disp(i)
            #disp(j)
    R = find_R(Table);
    simR = find_simR(R, Table);
    simR_Ax = find_SimR_Ax(simR, A.cT)
    simR_Bx = find_SimR_Ax(simR, B.cT)
    Score(i, j).lvalue = sum(simR_Ax) + sum(simR_Bx)
    Score(j, i).lvalue = sum(simR_Ax) + sum(simR_Bx)
    return(Score)

def find_R(Data=None):
    #UNTITLED2 Summary of this function goes here
    #   Detailed explanation goes here

    [m, n] = size(Data)
    R = zeros(n, m, m)
    for k in mslice[1:n]:
        R_k = find_Rk(k, Data)
        R(k, mslice[:], mslice[:]).lvalue = R_k
        end
        end
def find_Rk(k=None, Data=None):
    # This finction finds Lukasiewicz's opearator between two samples x_i and
    # x_j for the k^th gene
    # Data is the data matrix

    n = size(Data)
    for i in mslice[1:n]:
        for j in mslice[1:n]:
            R_k(i, j).lvalue = 1 - abs(Data(i, k) - Data(j, k))
            end
            end
            end
#@mfunction("simR_Ax")
def find_SimR_Ax(simR=None, A=None):
    # This function finds the 'C positive Region of D'
    # Here D = A
    r = size(A, 1)
    for i in mslice[1:r]:
        neg_R(i, mslice[:]).lvalue = 1 - simR(A(i), mslice[:])
        end
        s = size(neg_R, 2)
        B = mcat([mslice[1:s]])
        C = setdiff(B, A)
        s = size(C, 2)
        temp_R = zeros(r, s)
        for i in mslice[1:r]:
            for j in mslice[1:s]:
                temp_R(i, j).lvalue = neg_R(i, C(j))            # finds the matrix containing ~A set columns
                end
                end
                for i in mslice[1:r]:
                    simR_Ax(i).lvalue = min(temp_R(i, mslice[:]))
                    end
                    end
#@mfunction("simR")
def find_simR(R=None, Data=None):
    # R is the 3d matrix; R(k,i,j) = similarity between ith and jth patient for kth gene;

    [m, n] = size(Data)
    # n = number of genes

    for i in mslice[1:m]:
        for j in mslice[1:m]:
            simR(i, j).lvalue = min(R(mslice[:], i, j))
            end
            end
            end
