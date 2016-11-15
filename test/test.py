import numpy as np

a = np.array(
    [[8.1,2.1,3.2],
      [2.1,3.9,4.1],
      [2.1,4.1,1.1]]                   
    )

b = np.array(
[[8.1+2.1j,   6.2+12.1j,   6.2+12.1j],
[2.1+3.9j,  40.1+3.90j, 101.1+3.90j], 
[2.1+4.1j,   1.1+21.0j , 8.1+1.1j]] 
)

print np.linalg.det(a)
print np.linalg.det(b)
















































