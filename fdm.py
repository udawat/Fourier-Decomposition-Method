def fdm(X, fs, fc, data_type='columns', filter_type='dct', sort_fc='descend', remove_mean=False, plot_subbands=True):
# Take care of Vector inputs for einsum()
# ValueError: einstein sum subscripts string contains too many subscripts for operand 0
    if data_type == 'rows':
        # axis = 1
        X = np.transpose(X)
    else:
        axis = 0
    if remove_mean:
        X = X - np.mean(X, axis=0)
    N = X.shape[0]
    # print(N)
    # fc=np.asarray([0, fs/64, fs/32, fs/16, fs/8, fs/4, fs/2])
    fc = np.sort(fc)  # always sort fc in ascending order
    # print(fc)
    if fc[0] != 0:
        fc = np.hstack((0, fc))
    if fc[-1] != fs/2:
        fc = np.hstack((fc, fs/2))
    
    if sort_fc == 'descend':
        fc[::-1].sort()
    # print(fc)
    # Padding - Append data as symmetric rows (mirroring) at the beginning and end of the matrix
    # if filter_type != 'dct':
        # pad_rows = int(N * append_ratio);  # one-side length of appended rows
        # # X_padded = np.pad(X, pad_width=((pad_rows,pad_rows), (0,0)), mode="symmetric", reflect_type='even', )
        # X_padded = np.pad(X, pad_width=((pad_rows,pad_rows), (0,0)), mode="reflect", reflect_type='even')
        # appended_length = X_padded.shape[0]  # Increased number of rows before and after appending data, twice of 'pad_rows'

    # if filter_type == 'dct' or filter_type == 'dft':
    if filter_type == 'dct':
        dct_type = 2
        K = np.round(2 * N * fc/fs).astype(int) # Vector operation - convert fc to k (digital freq)
        no_of_subbands = K.shape[0] - 1         # For DCT implementation
        Hk = np.zeros((N, 1, no_of_subbands))   # Create 3-D Matrix with third dimension as Fc filters

    if filter_type == 'dft':
        append_ratio = 0.02  # 2 %
        pad_rows = int(N * append_ratio)  # one-side length of appended rows
        X_padded = np.pad(X, pad_width=((pad_rows,pad_rows), (0,0)), mode="symmetric", reflect_type='even')
        # X_padded = np.pad(X, pad_width=((pad_rows,pad_rows), (0,0)), mode="reflect", reflect_type='even')
        appended_length = X_padded.shape[0]  # Increased number of rows before and after appending data, twice of 'pad_rows'
        L = appended_length
        N_fft = 2 * np.ceil(L/2).astype(int)        # make it even
        K = np.round(N_fft * fc/fs).astype(int)     # Vector operation - convert fc to k (digital freq)
        # print(K)
        no_of_subbands = K.shape[0] - 1             # For FFT implementation
        Hk = np.zeros((N_fft, 1, no_of_subbands))   # Create 3-D Matrix with third dimension as Fc filters

    for i in range(no_of_subbands):
        if filter_type == 'dct':
            if sort_fc == 'ascend':
                Hk[K[i] : K[i+1], :, i] = 1  # 0:48 slice means 0 to 47 only
            if sort_fc == 'descend':
                Hk[K[i+1] : K[i], :, i] = 1
        if filter_type == 'dft':
            if sort_fc == 'ascend':
                Hk[K[i] : K[i+1], :, i] = 1
                Hk[N_fft-K[i+1] : N_fft-K[i], :, i] = 1
            if sort_fc == 'descend':
                Hk[K[i+1] : K[i], :, i] = 1
                Hk[N_fft-K[i] : N_fft-K[i+1], :, i] = 1
    
    if filter_type == 'dct':
        Xk = dct(X, type=dct_type, n=N, axis=axis, norm='ortho', overwrite_x=False, workers=None)
        # if X.shape[1] == 1:
        #     Yk = Xk * Hk
        # else:
        #     Yk = np.einsum('ij,ijk->ijk', Xk, Hk)
        Yk = np.einsum('ij,ijk->ijk', Xk, Hk)  # Python Broadcasting, Order of Multiplication of dimensions
        Y = idct(Yk, type=dct_type, n=N, axis=axis, norm='ortho', overwrite_x=False, workers=None)
    
    if filter_type == 'dft':
        Xk = 1/L * fft(X_padded, n=N_fft, axis=axis, norm=None)  # Divide by L because of scaling in FFT & IFFT formulae
        # if X.shape[1] == 1:
        #     Yk = Xk * Hk
        # else:
        #     Yk = np.einsum('ij,ijk->ijk', Xk, Hk)  # Python Broadcasting
        Yk = np.einsum('ij,ijk->ijk', Xk, Hk)  # Python Broadcasting
        Y = L * ifft(Yk, n=N_fft, axis=axis, norm=None, overwrite_x=False, workers=None)
        Y = np.real(Y)  # Ignore the imaginary part
        Y = Y[pad_rows+1 : appended_length-pad_rows+1, :, :]  # Remove padded data i.e., revert back to orginal size of data
        
    if X.shape[1] == 1:  # Output is returned as a 2D matrix if input to the function is a column vector or a row vector
        FIBFs = np.squeeze(Y)
    else:
        FIBFs = Y  # Return output in 3D Matrix with the third dimension corresponding to different subbands    

    if plot_subbands:
        # plt.show()
        t = np.arange(0, X.shape[0], 1) / fs
        fc = np.sort(fc)  # always sort fc in ascending order
        if fc[0] != 0:
            fc = [0, fc]
        if fc[-1] != fs/2:
            fc = [fc, fs/2]
        
        if sort_fc == 'descend':
            fc[::-1].sort()

        no_of_subbands = len(fc) - 1
        if no_of_subbands >= 6:
            m = np.floor(no_of_subbands/2) + 1
            n = 2  # Change to two-column tiled layout
        else:
            m = no_of_subbands + 1
            n = 1

        # First create a grid of plots
        # ax will be an array of two Axes objects
        # fig, ax = plt.subplots(m, n)
        # fig, ax = plt.subplots(m, n, figsize = (15, 10))
        # Call plot() method on the appropriate object
        # ax[0].plot(x, np.sin(x))
        # ax[1].plot(x, np.cos(x));
        if n == 2:
            plt.figure(figsize=(16, 16))
        else:
            plt.figure(figsize=(8, 16))
        # print(m)
        # print(n)
        
        for k in range(no_of_subbands+1):
            plt.subplot(m, n, k+1)
            if k == 0:
                plt.plot(t, X[:, 0])
                plt.title('Signal')
            else:
                if X.shape[1] == 1:  # Check if column vector or matrix
                    plt.plot(t, FIBFs[:, k-1])  # Plotting just one channel of the multi channel data
                else:
                    plt.plot(t, FIBFs[:, 0, k-1])  # Plotting just one channel of the multi channel data
                
                if sort_fc == 'ascend':
                    plt.title('FIBF {0}: Between {1} Hz to {2} Hz'.format(k, fc[k-1], fc[k]))
                elif sort_fc == 'descend':
                    plt.title('FIBF {0}: Between {1} Hz to {2} Hz'.format(k, fc[k], fc[k-1]))
        # plt.title('Fourier Decomposition Method')

    if data_type == 'rows':
        if X.shape[1] == 1:  # Output is returned as a 2D matrix if input to the function is a column vector or a row vector
            return FIBFs.T
        else:
            return FIBFs.transpose(1, 0, 2)  # The final 3D Matrix is a collection of multiple 2D matrices
    return FIBFs
