package matrix

import scala.collection.mutable.Buffer
import scala.math.*

/** A class of (mathematical) matrix objects of dimensions m*n containing Double values.
 * @param rows The number of rows (m) that the newly created matrix will have. Must be an integer greater than 0.
 * @param columns The number of columns (n) that the newly created matrix will have. Must be an integer greater than 0.
 * @param collection A collection of Double-type values to fill the matrix with, supplied in a single object.
 *                   If there are fewer objects than slots in the matrix (that is, rows * columns), the missing values will default to 0.0.
 *                   Conversely, if there are too many elements, only the first (rows * columns) will be taken into account.
 */
class Matrix(val rows: Int, val columns: Int, collection: scala.collection.Seq[Double]):

  private val slots = rows * columns
  private val itemCollection = collection.toVector.take(slots) ++ Vector.fill(slots - collection.size)(0.0)
  private val splicedCollection = itemCollection.grouped(columns).toVector

  //Helper method that returns a Vector of elements in specified row.
  //This is to allow more natural referral to rows in other methods, like in mathematics; that is, the first row is row 1, not 0.
  private def getRow(number: Int): Vector[Double] = this.splicedCollection(number - 1)

  //Helper method that performs an operation which takes an element from the first collection,
  //and an element with the same index from another collection, repeating this behaviour for all indices.
  private def operateOnMatching(collection1: Vector[Double], collection2: Vector[Double],
                                operation: (Double, Double) => Double): Vector[Double] =
    Vector.tabulate(collection1.size)(index => operation(collection1(index), collection2(index)))

  //Helper method that returns a Vector, created by performing a given operation on elements of two rows with matching indices.
  private def rowOperations(row1: Int, row2: Int, operation: (Double, Double) => Double): Vector[Double] =
    operateOnMatching(this.getRow(row1), this.getRow(row2), operation)

  //Helper method that returns a new Matrix object with contents of a specified row replaced with those of a given Vector.
  private def updateRow(rowToUpdate: Int, rowContents: Vector[Double]): Matrix =
    Matrix(this.rows, this.columns, this.splicedCollection.updated(rowToUpdate - 1, rowContents).flatten)

  //Helper method that returns a Vector, formed of elements of a given row multiplied by a given scalar value.
  private def multipliedRowVector(row: Int, scalar: Double): Vector[Double] =
    this.getRow(row).map(_ * scalar)

  /** Returns the number value slots in the matrix, that is: its row count multiplied by its column count.
   * @return An Int descrbing how many elements can this matrix store.
   */
  def slotCount: Int = this.slots

  /** Finds the n-th element in a matrix, where n is the index supplied as a parameter.
   * @param index A literal (that is, independent of the row and column count) index of the queried element.
   * @return Element from the matrix that corresponds to the provided index.
   */
  def element(index: Int): Double = this.itemCollection(index)

  /** Returns an element from the matrix, which can be found under the intersection of the specified row and column.
   * @param row The row in which the queried element is in.
   * @param column The column in which the queried element is in.
   * @return Element from the matrix that occupies the specified slot.
   */
  def element(row: Int, column: Int): Double = this.itemCollection((row - 1) * this.columns + column - 1)

  /** A way to access the private `val` containing all elements in the matrix in a one-dimensional vector form.
   * @return A vector containing all the values in this matrix.
   */
  def elementsVector: Vector[Double] = this.itemCollection

  /** A way to access the private `val` containing all elements in the matrix in a two-dimentional vector form.
   * @return A vector consisting of subvectors with all the values in this matrix, divided as per the matrix's rows.
   */
  def rowsVector: Vector[Vector[Double]] = this.splicedCollection

  /** A way to obtain all elements in the matrix in a two-dimentional vector form, grouped by columns they're stored in within the matrix.
   *
   * Same as calling rowsVector on a transposed matrix.
   * @see rowsVector
   * @return A vector consisting of subvectors with all the values in this matrix, divided as per the matrix's columns.
   */
  def columnsVector: Vector[Vector[Double]] = this.transpose.splicedCollection

  /** Describes the order (dimensions) of this matrix in a string form. Used by toString.
   * @see toString
   * @return A string in the format (mxn), where m and n are the matrix's row and column count, respectively.
   */
  def order: String = s"${this.rows}x${this.columns}"

  /** When called as a method of a Matrix object, performs the transposition operation on it - that is, reflects the matrix along its main diagonal.
   * @see [[https://en.wikipedia.org/wiki/Matrix_transposition]]
   * @example A 3x1 matrix becomes a 1x3 matrix with the same contents, but their row and column coordinates are swapped.
   * @return A new Matrix object, arranged so that the original matrix's contents are reflected along the main diagonal.
   */
  def transpose: Matrix =
    val outputCollection = for
      n <- 0 until this.columns
      m <- 0 until this.rows
    yield this.splicedCollection(m)(n)
    Matrix(rows = this.columns, columns = this.rows, outputCollection)

  /** Multiples all elements in a matrix by a given (scalar) value.
   * @see multiply(matrix: Matrix)
   * @param scalar The value by which all elements should be multipled.
   * @return A new Matrix object where all elements are scaled by the `scalar` factor.
   */
  def multiply(scalar: Double): Matrix =
    Matrix(this.rows, this.columns, this.itemCollection.map(_ * scalar))

  /** Performs the operation of matrix multiplication on two Matrix objects, with the first in order being the one upon this method is called.
   *
   * If the first matrix is of dimentsions m x n, the other matrix has to be of dimensions n x p.
   * @see [[https://en.wikipedia.org/wiki/Matrix_multiplication]]
   * @param matrix A matrix by which the first one should be multipled. Has to have as many rows as the first matrix has columns.
   * @return A new Matrix object of order m x p (same number of rows as the first matrix, and same number of columns as the second matrix).
   */
  def multiply(matrix: Matrix): Matrix =
    val outputCollection = ( for
      rowA <- 0 until this.rows
      rowB <- 0 until matrix.columns
    yield (0 until this.columns).foldLeft(0.0)( (sum: Double, item: Int) =>
          sum + this.splicedCollection(rowA)(item) * matrix.columnsVector(rowB)(item) )
      )
    Matrix(this.rows, matrix.columns, outputCollection)

  /** Shorthand notation identical to `multiply(scalar: Double)`.
   * @see multiply(scalar: Double)
   */
  def *(scalar: Double): Matrix = this.multiply(scalar)

  /** Shorthand notation identical to `multiply(matrix: Matrix)`.
   * @see multiply(matrix: Matrix)
   */
  def *(matrix: Matrix): Matrix = this.multiply(matrix)

  /** Performs the operation of matrix addition on two Matrix objects - that is, sums the elements from both matrices which have the same indices.
   *
   * The two matrices have to have the same dimensions.
   * @see [[subtract]]
   * @param matrix A matrix which should be added to the one upon which this method is called. Has to have the same dimensions as the first matrix.
   * @return A new Matrix object, with each element being a sum of the corresponding ones in the input matrices.
   */
  def add(matrix: Matrix): Matrix =
    Matrix(this.rows, this.columns, this.operateOnMatching(this.itemCollection, matrix.itemCollection, _ + _))

  /** Performs the operation of matrix subtraction on two Matrix objects -
   * that is, subtracts the elements from the second matrix from matching elements in the matrix upon which this method is called.
   *
   * The two matrices have to have the same dimensions.
   * @see [[add]]
   * @param matrix A matrix which should be subtracted from the one upon which this method is called. Has to have the same dimensions as the first matrix.
   * @return A new Matrix object, with each element being the difference of corresponding items in the input matrices.
   */
  def subtract(matrix: Matrix): Matrix =
    Matrix(this.rows, this.columns, this.operateOnMatching(this.itemCollection, matrix.itemCollection, _ - _))

  /** Multiplies elements with the same index found in two matrices by one another.
   *
   * The two matrices have to have the same dimensions.
   * @see [[divideMatching]]
   * @example If called upon a `[0.0, 1.0, 2.0]` matrix with `[3.0, 2.0, 1.0]` given as the parameter, it returns a `[0.0, 2.0, 2.0]` matrix.
   * @param matrix A matrix containing the elements by which those in the first matrix (one that this function is called upon) should be multiplied by.
   * @return A new Matrix object, with each element being the product of corresponding items in the input matrices.
   */
  def multiplyMatching(matrix: Matrix): Matrix =
    Matrix(this.rows, this.columns, this.operateOnMatching(this.itemCollection, matrix.itemCollection, _ * _))

  /** Divides the elements in the matrix upon which this method is called, by elements with a matching index found in the matrix supplied as a parameter.
   *
   * The two matrices have to have the same dimensions.
   * @see [[multiplyMatching]]
   * @example If called upon a `[0.0, 1.0, 2.0]` matrix with `[3.0, 2.0, 1.0]` given as the parameter, it returns a `[0.0, 0.5, 2.0]` matrix.
   * @param matrix A matrix containing the elements by which those in the first matrix (one that this function is called upon) should be divided by.
   * @return A new Matrix object, with each element being the quotient of corresponding items in the input matrices.
   */
  def divideMatching(matrix: Matrix): Matrix =
    Matrix(this.rows, this.columns, this.operateOnMatching(this.itemCollection, matrix.itemCollection, _ / _))

  /** Returns the remainders of dividing each element in the matrix that the method is called upon, by the elements with same index found in the matrix supplied as a parameter.
   *
   * The two matrices have to have the same dimensions.
   * @see [[divideMatching]]
   * @example If called upon a `[0.0, 1.0, 2.0]` matrix with `[3.0, 2.0, 1.0]` given as the parameter, it returns a `[0.0, 1.0, 0.0]` matrix.
   * @param matrix A matrix containing the elements by which those in the first matrix (one that this function is called upon) should be divided by.
   * @return A new Matrix object, with each element being the remainder of dividing items in the input matrices.
   */
  def remainderOfMatching(matrix: Matrix): Matrix =
    Matrix(this.rows, this.columns, this.operateOnMatching(this.itemCollection, matrix.itemCollection, _ % _))

  /** Shorthand notation identical to [[add]].
   * @see [[add]]
   */
  def +(matrix: Matrix): Matrix = this.add(matrix)

  /** Shorthand notation identical to [[subtract]].
   * @see [[subtract]]
   */
  def -(matrix: Matrix): Matrix = this.subtract(matrix)

  /** Performs the row operation of adding elements from one row to the elements in another row, replacing the other row in place.
   *
   * @see [[rowSubtract]]
   * @param addFromRow The row whose values should be added to the values found in the other row.
   * @param addAtRow The row which will be replaced with sums of the original values and the values in the first row (specified in the first parameter) with matching indices.
   * @return A new Matrix object with the row specified in the second parameter updated in place.
   */
  def rowAdd(addFromRow: Int, addAtRow: Int): Matrix =
    this.updateRow(addAtRow, this.rowOperations(addFromRow, addAtRow, _ + _))

   /** Performs the row operation of adding elements from one (scaled) row to the elements in another row, replacing the other row in place.
   *
   * @see [[rowSubtract]]
   * @param scalar A scale factor for the values found in the first row.
   * @param addFromRow The row whose scaled values should be added to the values found in the other row.
   * @param addAtRow The row which will be replaced with sums of the original values and the values in the first row (specified in the first parameter) with matching indices.
   * @return A new Matrix object with the row specified in the second parameter updated in place.
   */
  def rowAdd(scalar: Double, addFromRow: Int, addAtRow: Int): Matrix =
   this.updateRow(addAtRow, this.operateOnMatching(splicedCollection(addAtRow - 1), multipliedRowVector(addFromRow, scalar), _ + _))

  /** Performs the row operation of subtracting elements from one row from the elements in another row, replacing the other row in place.
   *
   * @see [[rowAdd]]
   * @param subtractAtRow The row which will be replaced with the difference of its original values and the values in the row specified in the second parameter with matching indices.
   * @param subtractWithRow The row whose values should be subtracted from the values found in the other specified row.
   * @return A new Matrix object with the row specified in the second parameter updated in place.
   */
  def rowSubtract(subtractAtRow: Int, subtractWithRow: Int): Matrix =
    this.updateRow(subtractAtRow, this.rowOperations(subtractAtRow, subtractWithRow, _ - _))


  /** Performs the row operation of subtracting elements from one (scaled) row from the elements in another row, replacing the other row in place.
   *
   * @see [[rowAdd]]
   * @param scalar A scale factor for the values found in the first row.
   * @param subtractAtRow The row which will be replaced with the difference of its original values and the values in the row specified in the second parameter with matching indices.
   * @param subtractWithRow The row whose scaled values should be subtracted from the values found in the other row.
   * @return A new Matrix object with the row specified in the second parameter updated in place.
   */
  def rowSubtract(scalar: Double, subtractAtRow: Int, subtractWithRow: Int): Matrix =
    this.updateRow(subtractAtRow, this.operateOnMatching(splicedCollection(subtractAtRow - 1), multipliedRowVector(subtractWithRow, scalar), _ - _))

  /** Updates the elements of a given row with the product of their multiplication by a specified value.
   *
   * @see [[rowDivide]]
   * @param row The row that should have its elements multiplied by the scalar.
   * @param scalar The value by which all elements in the specified row should be multipled.
   * @return A new Matrix object with the row specified in the first parameter scaled in place by the factor found in the second parameter.
   */
  def rowMultiply(row: Int, scalar: Double): Matrix =
    this.updateRow(row, multipliedRowVector(row, scalar))

  /** Updates the elements of a given row with the quotient of their division by a specified value.
   *
   * @see [[rowMultiply(]]
   * @param row The row that should have its elements divided by the scalar.
   * @param scalar The value by which all elements in the specified row should be divided.
   * @return A new Matrix object with the row specified in the first parameter divided in place by the value found in the second parameter.
   */
  def rowDivide(row: Int, scalar: Double): Matrix = rowMultiply(row, 1 / scalar)

  /** Switches two specified rows with one another.
   *
   * @param row1 The first row to be swapped.
   * @param row2 The second row to be swapped.
   * @return A new Matrix object with the two specified rows swapped with one another.
   */
  def rowSwitch(row1: Int, row2: Int): Matrix =
    val oldRow1Vector = this.getRow(row1)
    this.updateRow(row1, this.getRow(row2)).updateRow(row2, oldRow1Vector)

  /** Creates a new, smaller matrix that is an intersection of select rows and columns found in the matrix upon which this function is called.
   *
   * Mathematically, all values in the output matrix should be adjacent to one another in the original matrix, but this is not checked to allow possible other, unorthodox uses.
   * @param rowsToTake        A Vector object listing the rows which have at least one element that should be included in the output matrix.
   * @param columnsToSubtract A Vector object listing the columns which should NOT be included in the output matrix.
   * @return A new Matrix object that is formed of the values found in the given rows belonging to the original matrix, with the columns specified in the second parameter omitted.
   */
  def submatrix(rowsToTake: scala.collection.Seq[Int], columnsToSubtract: scala.collection.Seq[Int]): Matrix =
    val newRows = rowsToTake.map(this.getRow(_).map(Option(_)).toBuffer)
    for (row <- newRows.indices; column <- columnsToSubtract) do
      newRows(row)(column - 1) = None
    val outputRows = newRows.flatMap(_.filter(_.isDefined).flatten)
    Matrix(rowsToTake.length, (this.columns - columnsToSubtract.length), outputRows)

  /** Checks if the matrix upon which this method is called is in row echelon form.
   *
   * @see [[https://en.wikipedia.org/wiki/Row_echelon_form]]
   * @return A Boolean value based on whether the matrix passed this check.
   */
  def isInRowEchelonForm: Boolean =
    var successFlag = true
    for column <- 0 until (this.columns - 1) do
      successFlag = successFlag && (columnsVector(column).drop(column + 1).sum == 0)
    for row <- 0 until this.rows do
      if this.splicedCollection(row).sum == 0 then
        successFlag = successFlag && (row == this.rows - 1 || this.splicedCollection.map(_.sum).drop(row + 1).sum == 0)
    successFlag

  /** Obtains the trace of the (square) matrix upon which this method is called. This is the same as the determinant of a matrix in row echelon form.
   *
   * Requires the matrix upon which this method is called to be square; that is, have the same number of rows as columns.
   * @see [[determinant]]
   * @see [[https://en.wikipedia.org/wiki/Matrix_(mathematics)#Trace]]
   * @see [[https://en.wikipedia.org/wiki/Row_echelon_form]]
   * @return The trace value wrapped in an Option, if the matrix is square; `None` otherwise.
   */
  def trace: Option[Double] =
    var output: Option[Double] = None
    var rawOutput = 0.0

    if this.rows == this.columns then
      for index <- 0 until this.rows do
        rawOutput += this.splicedCollection(index)(index)
      output = Some(rawOutput * this.echelon(1))
    output
  end trace

  /** Creates a new matrix by putting the matrix upon which this method is called to the left of the matrix supplied as argument.
   * That is, each row of the first matrix has the matching rows of the second matrix appended.
   *
   * @param another The matrix whose rows' contents should be appended to this matrix's rows.
   * @return A new matrix of order m x n+p, where:
   *         m is the number of rows (same for both original matrices),
   *         n is the number of columns of the first matrix, and
   *         p is the number of columns of the second matrix.
   * @throws IllegalArgumentException if there is a mismatch between the number of rows in the original matrices.
   * @see [[rightOf]]
   *      [[above]]
   *      [[below]]
   */
  infix def leftOf(another: Matrix): Matrix =
    if this.rows != another.rows then throw IllegalArgumentException("The row counts of the two matrices supplied don't match.")
    val newSpliced = this.splicedCollection.indices.map( index => this.splicedCollection(index) ++ another.splicedCollection(index) )
    Matrix(this.rows, (this.columns + another.columns), newSpliced.flatten)

  /** Creates a new matrix by putting the matrix upon which this method is called to the right of the matrix supplied as argument.
   * That is, each row of the second matrix has the matching rows of the first matrix appended.
   *
   * @param another The matrix whose rows should have this matrix's rows' contents appended to them.
   * @return A new matrix of order m x n+p, where:
   *         m is the number of rows (same for both original matrices),
   *         n is the number of columns of the first matrix, and
   *         p is the number of columns of the second matrix.
   * @throws IllegalArgumentException if there is a mismatch between the number of rows in the original matrices.
   * @see [[leftOf]]
   *      [[above]]
   *      [[below]]
   */
  infix def rightOf(another: Matrix): Matrix =
    another.leftOf(this)

  /** Creates a new matrix by putting the matrix upon which this method is called above the matrix supplied as argument.
   * That is, each column of this matrix has the matching columns of the second matrix appended.
   *
   * @param another The matrix whose columns' contents should be appended to this matrix's columns.
   * @return A new matrix of order m x n+p, where:
   *         m is the number of rows (same for both original matrices),
   *         n is the number of columns of the first matrix, and
   *         p is the number of columns of the second matrix.
   * @throws IllegalArgumentException if there is a mismatch between the number of rows in the original matrices.
   * @see [[leftOf]]
   *      [[rightOf]]
   *      [[below]]
   */
  infix def above(another: Matrix): Matrix =
    this.transpose.leftOf(another.transpose).transpose

  /** Creates a new matrix by putting the matrix upon which this method is called below the matrix supplied as argument.
   * That is, each column of the other matrix has the matching columns of this matrix appended.
   *
   * @param another The matrix whose columns should have this matrix's columns' contents appended to them.
   * @return A new matrix of order m x n+p, where:
   *         m is the number of rows (same for both original matrices),
   *         n is the number of columns of the first matrix, and
   *         p is the number of columns of the second matrix.
   * @throws IllegalArgumentException if there is a mismatch between the number of rows in the original matrices.
   * @see [[leftOf]]
   *      [[rightOf]]
   *      [[above]]
   */
  infix def below(another: Matrix): Matrix =
    another.transpose.leftOf(this.transpose).transpose

  /** Obtains the determinant of the (square) matrix upon which this method is called. This is the same as the trace of a matrix in row echelon form.
   *
   * Requires the matrix upon which this method is called to be square; that is, have the same number of rows as columns.
   * @see [[trace]]
   * @see [[https://en.wikipedia.org/wiki/Determinant]]
   * @see [[https://en.wikipedia.org/wiki/Row_echelon_form]]
   * @return The trace value wrapped in an Option, if the matrix is square; `None` otherwise.
   */
  def determinant: Option[Double] =
    var output: Option[Double] = None
    var rawOutput = 0.0

    if this.rows == this.columns then
      for index <- 0 until this.rows do
        rawOutput *= this.echelon(0).splicedCollection(index)(index)
      output = Some(rawOutput * this.echelon(1))

    output
  end determinant

  /** Shorthand notation identical to [[determinant]].
   * @see [[determinant]]
   */
  def det: Option[Double] = this.determinant

  /** Work in progress - not implemented fully.
   * Does nothing in its current form except returning: the same matrix as the one upon which this method is called, if [[isInRowEchelonForm]] is true, and the number 1.
   *
   * To get the correct determinant of the output matrix, multiply the normally obtained result by the second element of the Vector returned by this function.
   *
   * This is because each row swapping operation impacts the determinant of a matrix by changing its sign. Thus, the second value returned by this function keeps track of
   * the number of such operations that were performed by the function whilst getting the original matrix to the row echelon form.
   * @see [[isInRowEchelonForm]]
   * @see [[https://en.wikipedia.org/wiki/Row_echelon_form]]
   * @see [[https://sites.ualberta.ca/~jsylvest/books/DLA/section-det-by-row-red-concepts.html]]
   * @return A new Matrix object in row echelon form, and an Int either equal to 1 or -1.
   */
  def echelon: (Matrix, Int) = //TODO: fix / rewrite
    if isInRowEchelonForm then
      (this, 1)
    else
      ???

//  loop N through columns
//  choose element N,N as pivot
//  loop E through rows starting at N + 1
//    if element E,N != 0 then
//      subtract (E,N) / (N,N) * (row N) from row E replacing row E
//

  /** Work in progress - not implemented.
   *
   * @see [[echelon]]
   * @return
   */
  def reducedEchelon: (Matrix, Int) =
    var operatedOn = (this, 1)
    if !isInRowEchelonForm then
      operatedOn = this.echelon
    else ???
    operatedOn

  /** Work in progress - not fully implemented.
   *
   * @return A new Matrix object wraped in an Option, which is the inverse of the original matrix (if the original matrix has an inverse).
   *         Otherwise, `None`.
   */
  def inverse: Option[Matrix] =
    if this.determinant.getOrElse(0) == 0 then None
    else
      var output = this.leftOf(IdentityMatrix(this.rows)).reducedEchelon(0)
      output = output.submatrix(Vector.tabulate(this.rows)(_ + 1), Vector.tabulate(this.columns)(_ + 1))
      Some(output)


  //TODO: Complete methods: inverse, echelon, reducedEchelon

  /** Returns a String that states the [[order]] (dimensions) of the matrix at hand and visualises it, with each new line being subsequent rows of the matrix.
   *
   * The width of each column scales dynamically depending on the length of its contents, so that only the elements in same columns are directly underneath each other.
   * Two spaces separate each column, and one space separates the first and last columns from the respective matrix boundaries marked with `[` and `]` symbols.
   *
   * @example ```text
   *          Matrix (3x2)
   *
   *          [ 2.0  3.14159  0.0 ]
   *          [ 1.0  3.0      7.0 ] ```
   *
   * @see [[order]]
   * @return A String object that is a visual description of the matrix.
   */
  override def toString: String =
    var output = s"Matrix (${this.order})"
    val maxLengthInEachColumn = Buffer[Int]()
    for column <- columnsVector do
      maxLengthInEachColumn += column.map(_.toString.length).max

    for row <- 0 until this.rows do
      output += "\n[ "
      for column <- 0 until this.columns do
        val element = this.itemCollection(row * this.columns + column).toString
        output += element
        output += " " * (maxLengthInEachColumn(column) - element.length) + "  "
      output += "\b]"

    output
  end toString


private def getIdentityMatrixContents(dimensions: Int): Vector[Double] =
  val repeating = Vector.tabulate(dimensions + 1)(index => if index % (dimensions + 1) == 0 then 1.0 else 0.0 )
  (1 until dimensions).flatMap(_ => repeating).toVector :+ 1.0

class IdentityMatrix(dimensions: Int) extends Matrix(dimensions, dimensions, getIdentityMatrixContents(dimensions))