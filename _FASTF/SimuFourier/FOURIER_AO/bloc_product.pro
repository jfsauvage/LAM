FUNCTION bloc_product, array1, array2, blocsize

; performs a matrix-matrix product, bloc by bloc. The bloc is assumed to be
; squared, but the matrices may be rectangular.

; array1 : first matrix 
; array2 : second matrix 
; blocsize : ARRAY1 and ARRAY2 are assumed to be of shape
; ARRAY1 : size1x x blocsize  X  size1y x blocsize
; ARRAY1 : size2x x blocsize  X  size2y x blocsize

size1x = (size(array1))[1] / blocsize
size1y = (size(array1))[2] / blocsize
size2x = (size(array2))[1] / blocsize
size2y = (size(array2))[2] / blocsize

IF size1x NE size2y THEN print, 'Size error for BLOCK_PRODUCT'

array_out = complexarr(size2x * blocsize, size1y * blocsize)

FOR i = 0, size2x-1 DO $
   FOR j = 0, size1y-1 DO $
       FOR k = 0, size1x-1 DO $
           array_out(i*blocsize:(i+1)*blocsize-1, j*blocsize:(j+1)*blocsize-1) += $
    array1(k*blocsize:(k+1)*blocsize-1, j*blocsize:(j+1)*blocsize-1) * $
    array2(i*blocsize:(i+1)*blocsize-1,k*blocsize:(k+1)*blocsize-1 )

return, array_out

end
