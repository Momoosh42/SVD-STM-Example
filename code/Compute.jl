function make_binary_stoichiometric_matrix(stoichiometric_matrix)

    # initialize -
    (number_of_rows, number_of_cols) = size(stoichiometric_matrix)
    binary_stoichiometric_matrix = zeros(number_of_rows, number_of_cols)
    for row_index = 1:number_of_rows
        for col_index = 1:number_of_cols

            # grab value -
            old_value = stoichiometric_matrix[row_index, col_index]
            if (old_value !=0.0)
                binary_stoichiometric_matrix[row_index, col_index] = 1.0
            end
        end
    end

    # return -
    return binary_stoichiometric_matrix
end

function fractional_reconstruction_of_stoichiometric_matrix(stoichiometric_matrix,U,S,V)

    # what is the rank of S?
    rank_value = rank(diagm(S))
    (number_of_rows, number_of_cols) = size(stoichiometric_matrix)
    block_array = Array{Array{Float64,2},1}()

    fraction_block = zeros(number_of_rows, number_of_cols)
    for rank_index = 1:rank_value

        fraction_block = fraction_block + S[rank_index]*(U[:,rank_index]*transpose(V[:,rank_index]))
        push!(block_array,fraction_block)
    end

    return block_array
end
