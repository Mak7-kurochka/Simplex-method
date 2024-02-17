import numpy as np

class LPP:
    def __init__(self,direction):
        self.direction = direction
    def find_change(self,cj_zj):
        maxx_val = 0
        for i in range(len(cj_zj)):
            if cj_zj[i] > maxx_val:
                maxx_val = cj_zj[i]
                maxx_num = i
        return maxx_num
    def change_cb(self,column,matrix,b,c_new,c_b,zm_baz):
        minn_val = max(b)
        for i in range(len(b)):
            if matrix[i][column] > 0 and b[i]/matrix[i][column] < minn_val:
                minn_val = b[i]/matrix[i][column]
                minn = i
        c_b[minn] = c_new[column]
        zm_baz[minn] = column
        return [c_b,zm_baz]
    def built_first_table(self,c,restr):
        n = len(c)
        m = len(restr)
        A = np.empty((m,n))

        for i in range(m):
            A[i] = restr[i][0]

        I = np.identity(m)

        b = np.empty(m)
        for i in range(m):
           b[i] = restr[i][1]

        c_new = np.append(c,np.zeros(m))
        c_b = np.zeros(m)
        zm_baz = np.arange(len(c),len(c_new))
        z_j = np.zeros(n+m)
        cj_zj = c_new - z_j
        return self.go(A,I,b,c_b,c_new,zm_baz, cj_zj)

    def sympleks(self,c,restr):
        if self.direction == 'max':
            return self.symplex_max(c,restr)
        else:
            return 'Unfortunately this direction of optimization doesnt work yet'
    
    def biult_B(self,A,I,c_b,zm_baz):
        n = A.shape[1]
        m = len(c_b)
        B = np.zeros((m,m))
        matrix = np.concatenate((A,I.T),axis=1)
        i = 0
        while i < B.shape[0]:
            B[:, i] = matrix[:, zm_baz[i]]
            i += 1
        return B

    def go(self,A,I,b,c_b,cj,zm_baz,cj_zj):
        column = self.find_change(cj_zj)
        if column < A.shape[1]:
            res = self.change_cb(column,A,b,cj,c_b,zm_baz)
            c_b = res[0]
            zm_baz = res[1]
        else:
            res = self.change_cb(column,A,b,cj,c_b,zm_baz)
            c_b = res[0]
            zm_baz = res[1]
        B = self.biult_B(A,I,c_b,zm_baz)
        Binv = np.linalg.inv(B)
        Binv_A = np.matmul(Binv, A)
        Binv_b = np.matmul(Binv,b)
        cb_Binv_A = np.matmul(c_b,Binv_A)
        cb_Binv = np.matmul(c_b,Binv)
        cb_Binv_b = np.matmul(c_b, Binv_b)
        row = np.concatenate((cb_Binv_A,cb_Binv))
        cj_zj = cj - row
        if False in (cj_zj<=0):
            return self.go(A,I,b,c_b,cj,zm_baz,cj_zj)
        else:
            return f'FC:{cb_Binv_b}, x:{zm_baz}, number:{Binv_b}'

    def symplex_max(self,c, restr):
        n = len(c)
        m = len(restr)
        for i in range(m):
            if len(restr[i][0]) < n:
                _ = [0 for x in range(m-len(restr[i][0])-1)]
                restr[i][0].extend(_)
        var = np.empty(n)
#        restr_new = list(restr)
#        for i in range(m):
#            free_val = np.zeros(m)
#            free_val[i] = 1
#            restr_new[i][0].extend(free_val)
        limit = np.empty(m)
        for i in range(m):
           limit[i] = (restr[i][0]@var <= restr[i][1])
        return self.built_first_table(c,restr)

problem = LPP('max')
print(problem.sympleks([30,20], [[[2,1],1000], [[3,3],2400],[[1.5], 600]]))