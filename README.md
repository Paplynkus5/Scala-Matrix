A utility for matrix calculations in Scala. 
Documentation is available under `doc\matrix.html` as a Scaladoc.

## Features implemented and left to do:

| Feature | Implementation |
| --- | ----------- |
| Matrix object type | `Matrix(val rows: Int, val columns: Int, collection: Seq[Double])` |
| Describing a matrix | `toString`, `order` |
| Matrix addition, subtraction | `add`, `subtract`, `+`, `-` |
| Multiplying matrices | `multiply(matrix: Matrix)`, `*(matrix: Matrix)` |
| Row operations | `rowAdd`, `rowSubtract`, `rowMultiply`, `rowDivide`, `rowSwitch` |
| Matrix transposition | `transpose` |
| Trace calculation | `trace` |
| Methods returning elements stored as collections | `rowsVector`, `columnsVector`, `elementsVector` |
| Methods returning a specified element | `element(index: Int)`, `element(row: Int, column: Int)` |
| Creating submatrices | `submatrix` |
| Operations on a matrix's elements with another matrix's matching elements as additional parameters | `divideMatching`, `multiplyMatching`, `remainderOfMatching`, `add`, `subtract`, `+`, `-` |
| Scaling | `multiply(scalar: Double)`, `*(scalar: Double)` |
| Obtaining the determinant of a (square) matrix | `determinant`, `det` (both partially implemented) |
| Inverse of a matrix | `inverse` (not yet implemented) |
| Obtaining the row echelon form of a matrix | `isInRowEchelonForm`, `echelon` (not yet implemented) |
| Reduced row echelon form | `reducedEchelon` (not yet implemented) |
