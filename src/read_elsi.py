import struct
import numpy as np

# Read ELSI binary matrix, convert to dense array
def read_elsi_to_den(filename):
    mat = open(filename,"rb")
    data = mat.read()
    mat.close()

    # Get header
    start = 0
    end = 64
    header = struct.unpack("i"*16,data[start:end])

    # Matrix size
    dim = header[3]

    # Number of non-zero elements
    nnz = header[5]

    # Get column pointer
    start = end
    end += dim*4
    col_ptr = struct.unpack("i"*dim,data[start:end])
    col_ptr += (nnz+1,)
    col_ptr = np.array(col_ptr)

    # Get row index
    start = end
    end += nnz*4
    row_idx = struct.unpack("i"*nnz,data[start:end])
    row_idx = np.array(row_idx)

    # Get non-zero value
    start = end

    if header[2] == 0:
        # Real case
        end += nnz*8
        nnz_val = struct.unpack("d"*nnz,data[start:end])
    else:
        # Complex case
        end += nnz*16
        nnz_val = struct.unpack("d"*nnz*2,data[start:end])
        nnz_val_real = np.array(nnz_val[0::2])
        nnz_val_imag = np.array(nnz_val[1::2])
        nnz_val = nnz_val_real + 1j*nnz_val_imag

    nnz_val = np.array(nnz_val)
    i_val = 0

    if header[2] == 0:
        # Real case
        den = np.zeros((dim,dim))

        for i_col in range(dim):
            this_nnz = col_ptr[i_col+1]-col_ptr[i_col]

            for j_val in range(i_val,i_val+this_nnz):
                den[i_col,row_idx[j_val]-1] = float(nnz_val[j_val])

            i_val += this_nnz
    else:
        # Complex case
        den = np.zeros((dim,dim),dtype=complex)

        for i_col in range(dim):
            this_nnz = col_ptr[i_col+1]-col_ptr[i_col]

            for j_val in range(i_val,i_val+this_nnz):
                den[i_col,row_idx[j_val]-1] = complex(nnz_val[j_val])

            i_val += this_nnz

    return den
