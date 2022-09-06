CC=gcc
CFLAGS= -g
ifeq ($(shell uname), Linux)
	CHKFLAGS=-pthread -lcheck_pic -pthread -lrt -lm -lsubunit
else
	CHKFLAGS=-lcheck
endif

all: clean s21_matrix.a test linter

clean:
	@find . -type f -name "*.a" -delete
	@find . -type f -name "*.o" -delete
	@find . -type f -name "*.out" -delete
	@rm -rf test *.html *.gcda *.gcno gcov_report *.css
	@rm -rf test test.c

s21_matrix.o: s21_matrix.c
	@gcc -o s21_matrix.o $^ -c

s21_matrix.a: s21_matrix.o
	@ar crs $@ $^

test:
	@echo 'test'
	$(CC) $(CFLAGS) s21_matrix.c -c -o s21_matrix_test.o
	$(CC) $(CFLAGS) tests.c -c -o test_test.o
	$(CC) test_test.o s21_matrix_test.o -o $@ $(CHKFLAGS)
	./$@

gcov_report:
	@echo 'gcov_report'
	$(CC) $(CFLAGS) --coverage s21_matrix.c -c -o s21_matrix_coverage.o
	$(CC) $(CFLAGS) tests.c -c -o test_coverage.o
	$(CC) s21_matrix_coverage.o test_coverage.o -o $@ $(CHKFLAGS) -lgcov
	./$@
	gcovr --html --html-details -o report.html
	open report.html


linter:
	@cp ../materials/linters/CPPLINT.cfg ./
	@find . -name "*.c" -exec python3 ../materials/linters/cpplint.py {} \;
	@find . -name "*.h" -exec python3 ../materials/linters/cpplint.py {} \;
	@rm -f CPPLINT.cfg
