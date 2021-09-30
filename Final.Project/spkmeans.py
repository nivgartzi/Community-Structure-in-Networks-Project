import numpy as np
import pandas as pd
import math
import argparse
import sys
import csv
import mykmeanssp as mk


"""HW2"""
def receiving_user_inputs():
    reader = argparse.ArgumentParser()
    reader.add_argument("k", type=int)
    reader.add_argument("goal", type=str)
    reader.add_argument("file_name", type=str)
    arguments = reader.parse_args()
    return arguments


def distance_calculator(v1, v2):
    if (len(v1) != len(v2)):
        print("Invalid Input!")
    distance = 0
    for c in range(len(v1)):
        coordinates_d = (v1[c] - v2[c]) ** 2
        distance += coordinates_d
    return distance


def probability_calculator(i, D_lst, D_sum):
    prob = float ( (D_lst[i] / D_sum) )
    return prob

def init_input(argument):
    input_df = pd.read_csv(argument, header=None)
    return input_df


#prints the output for python flow
def print_required(first_k_centroids_indices,FINAL_centroids_in_numpy_converted_to_list):  
    for i in range(len( first_k_centroids_indices)):
        if (i != (len( first_k_centroids_indices) - 1)):
            print(  first_k_centroids_indices[i] , end=","  )
        else:
            print( first_k_centroids_indices[i])

    n=len(FINAL_centroids_in_numpy_converted_to_list)
    m=len(FINAL_centroids_in_numpy_converted_to_list[0])
    for i  in range(n-1):
        for j in range(m):
            if(j!=m-1):
             print(format(FINAL_centroids_in_numpy_converted_to_list[i][j],".4f"), end=",")
            else:
                print(format(FINAL_centroids_in_numpy_converted_to_list[i][j],".4f"), end="")
        print()   
    for j in range(m):
        if(j!=m-1):
         print(format(FINAL_centroids_in_numpy_converted_to_list[n-1][j],".4f"), end=",")
        else:
            print(format(FINAL_centroids_in_numpy_converted_to_list[n-1][j],".4f"), end="" )



def find_first_indices_and_corresponding_vectors(k, data_array, n, d):
    np.random.seed(0)
    first_centroid_index = np.random.choice(n)  #choose row number randomally
    initial_centroids_indices = [-1 for i in range(k+1)]
    initial_centroids_indices[1] = first_centroid_index
    D_lst = [0 for i in range(n)]
    P_lst = [0 for i in range(n)]
    z = 1

    first_k_centroids = np.empty((0, d))
    first_k_centroids = np.append(first_k_centroids, [data_array[first_centroid_index]], axis=0)

    while ( z < k ):
        D_sum = 0
        for i in range(n):   #for each xi
            min_distance = math.inf
            for j in range(1,z+1):
                v1 = data_array[i]
                indice_of_chosen_vector_in_j_th_folder = initial_centroids_indices[j]
                v2 = data_array[indice_of_chosen_vector_in_j_th_folder]
                curr_distance = distance_calculator(v1,v2)
                if(curr_distance<=min_distance):
                    min_distance=curr_distance
            D_lst[i] = min_distance
        for i in range(n):   #calculate D_sum
            D_sum += D_lst[i]
        for i in range(n):  #find P(Xi) for each Xi and update P_lst
            P_lst[i] = probability_calculator(i,D_lst,D_sum)
        z += 1
        z_centroid_index = np.random.choice(n,p=P_lst)
        initial_centroids_indices[z] = z_centroid_index
        first_k_centroids = np.append(first_k_centroids, [data_array[z_centroid_index]], axis=0)
    return first_k_centroids, initial_centroids_indices[1:len(initial_centroids_indices)]
""" """

def main_funct():
    #set k, goal and file content
    arguments = receiving_user_inputs()
    external_k = arguments.k
    goal = arguments.goal
    file_content_1_pandas = init_input(arguments.file_name)
    file_content_1_numpy = file_content_1_pandas.to_numpy()
    n = file_content_1_numpy.shape[0]
    d_of_file_content_1 = file_content_1_numpy.shape[1]
    file_content_1 = file_content_1_numpy.tolist()
    

    if(goal=="spk"):
        k = mk.getk( file_content_1 ,  n,  d_of_file_content_1 , external_k)
        if (k == -1):
            print("An Error Has Occured")

        file_content_2 = mk.get_file_content_2 (file_content_1 ,  n,  d_of_file_content_1 , external_k, k)

        #at file_content_2 the dimension is k
        d=k
        first_k_centroids, first_k_centroids_indices =  find_first_indices_and_corresponding_vectors(k, file_content_2, n, d)

        first_k_centroids = first_k_centroids.tolist()

        the_final_centroids = mk.get_final_clusters(first_k_centroids, file_content_2, 300, n, d, k)

        the_final_centroids_in_numpy = np.array(the_final_centroids)

        print_required(first_k_centroids_indices, the_final_centroids_in_numpy.tolist() )

    if(goal=="wam"):
        id = mk.other_cases(file_content_1,n, d_of_file_content_1,2)
    if(goal=="ddg"):
        id = mk.other_cases(file_content_1,n, d_of_file_content_1,3)
    if(goal=="lnorm"):
        id = mk.other_cases(file_content_1,n, d_of_file_content_1,4)
    if(goal=="jacobi"):
        id = mk.other_cases(file_content_1,n, d_of_file_content_1,5)

main_funct()













