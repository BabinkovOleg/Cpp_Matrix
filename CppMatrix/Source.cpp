#include <iostream>
#include <time.h>
#include <algorithm>

#define RAND_LIMIT 10
#define EPSILON 0.0000001f
#define ERR_ALLOC "Allocation failed"
#define ERR_OPERATOR "Invalid matrix dimensions"
#define ERR_COORDINATES "Invalid element coordinates"
#define ERR_DIV_BY_NULL "Dividing by zero" 

struct Matrix {
private:
	int rows;
	int columns;
	float** elements;
public:
	void MatrixAlloc() {
		this->elements = new float* [this->rows];
		if (!this->elements) {
			std::cout << ERR_ALLOC << std::endl;
			delete[] this->elements;
			exit(1);
		}

		for (int i = 0; i < this->rows; ++i) {
			this->elements[i] = new float[this->columns];
			if (!this->elements[i]) {
				for(int j = 0; j < i; ++j){
					delete[] this->elements[j];
				}
				delete[] this->elements;
				std::cout << ERR_ALLOC << std::endl;
				exit(1);
			}
		}
	}

	void input(int _rows, int _columns) {
		this->rows = _rows;
		this->columns = _columns;
		this->MatrixAlloc();

		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				std::cin >> this->elements[i][j];
			}
		}
	}

	void generateRandom(int _rows, int _columns) {
		this->rows = _rows;
		this->columns = _columns;
		this->MatrixAlloc();

		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				this->elements[i][j] = rand() % RAND_LIMIT;
			}
		}
	}

	void output() {
		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				std::cout << this->elements[i][j] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
	void plus(Matrix m1, Matrix m2) {
		if (m1.rows != m2.rows || m1.columns != m2.columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		this->rows = m1.rows;
		this->columns = m1.columns;
		this->MatrixAlloc();

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] = m1.elements[i][j] + m2.elements[i][j];
			}
		}
	}

	void minus(Matrix m1, Matrix m2) {
		if (m1.rows != m2.rows || m1.columns != m2.columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		this->rows = m1.rows;
		this->columns = m1.columns;
		this->MatrixAlloc();

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] = m1.elements[i][j] - m2.elements[i][j];
			}
		}
	}

	void multiplication(float number, Matrix m) {
		this->rows = m.rows;
		this->columns = m.columns;

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] = m.elements[i][j] * number;
			}
		}
	}

	void multiplication(Matrix m1, Matrix m2) {
		if (m1.columns != m2.rows) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		this->rows = m1.rows;
		this->columns = m2.columns;
		this->MatrixAlloc();

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] = 0;
				for (int k = 0; k < m1.columns; ++k) {
					this->elements[i][j] += m1.elements[i][k] * m2.elements[k][j];
				}
			}
		}
	}

	int getElement(int i, int j) {
		if (i > this->rows || j > this->columns || i < 0 || j < 0) {
			std::cout << ERR_COORDINATES << std::endl;
			exit(3);
		}

		return this->elements[i][j];
	}

	void swapRows(int i1, int i2) {
		for (int j = 0; j < this->columns; ++j) {
			std::swap(this->elements[i1][j], this->elements[i2][j]);
		}
	}

	float determinant() {
		if (this->rows != this->columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		Matrix m;
		m.rows = this->rows;
		m.columns = this->columns;
		m.MatrixAlloc();

		float result = 1;
		
		for (int i = 0; i < m.rows; ++i) {
			for (int j = 0; j < m.columns; ++j) {
				m.elements[i][j] = this->elements[i][j];
			}
		}

		for (int i = 0; i < m.columns - 1; ++i) {
			for (int j = i + 1; j < m.rows; ++j) {
				if (m.elements[i][i] == 0) {
					int l = i + 1;
					while (m.elements[l][i] == 0 && l < m.rows - 1) {
						++l;
					}
					if (m.elements[l][i] == 0) {
						return 0;
					}
					m.swapRows(i, l);
					result *= -1;
				}
				float coefficient = m.elements[j][i] / m.elements[i][i];
				for (int k = i; k < m.columns; ++k) {
					m.elements[j][k] -= coefficient * m.elements[i][k];
				}
			}
		}
		for (int i = 0; i < m.rows; ++i) {
			result *= m.elements[i][i];
		}
		return result;
	}

	void transpose(Matrix m) {
		this->rows = m.columns;
		this->columns = m.rows;
		this->MatrixAlloc();

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] = m.elements[j][i];
			}
		}
	}

	float minor(int i, int j) {
		if (this->rows != this->columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		Matrix m;
		m.rows = this->rows - 1;
		m.columns = this->columns - 1;
		m.MatrixAlloc();

		for (int k = 0; k < m.rows; ++k) {
			for (int l = 0; l < m.columns; ++l) {
				if (k >= i) {
					if (l >= j) {
						m.elements[k][l] = this->elements[k + 1][l + 1];
					}
					else {
						m.elements[k][l] = this->elements[k + 1][l];
					}
				}
				else {
					if (l >= j) {
						m.elements[k][l] = this->elements[k][l + 1];
					}
					else {
						m.elements[k][l] = this->elements[k][l];
					}
				}
			}
		}
		return m.determinant();
	}

	void invert(Matrix m) {
		if (m.rows != m.columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}
		if (abs(m.determinant() - 0.0f) < EPSILON) {
			std::cout << ERR_DIV_BY_NULL << std::endl;
			exit(4);
		}

		this->rows = m.rows;
		this->columns = m.columns;
		this->MatrixAlloc();

		Matrix tmp;
		tmp.adjugate(m);
		tmp.transpose(tmp);
		this->multiplication(1 / m.determinant(), tmp);
	}

	float cofactor(int i, int j) {
		return ((i + j) % 2 == 0 ? 1 : -1) * this->minor(i, j);
	}

	void adjugate(Matrix m) {
		if (m.rows != m.columns) {
			std::cout << ERR_OPERATOR << std::endl;
			exit(2);
		}

		this->rows = m.rows;
		this->columns = m.columns;
		this->MatrixAlloc();

		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->columns; ++j) {
				this->elements[i][j] =m.cofactor(i, j);
			}
		}
	}
};

int main() {
	srand(time(NULL));
	Matrix m, mInverted;
	m.generateRandom(3, 3);
	m.output();
	mInverted.invert(m);

	mInverted.output();
	return 0;
}