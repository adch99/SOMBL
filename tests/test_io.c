#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../extern/unity/unity.h"
#include "../src/io/io.h"
#include "../src/utils/utils.h"
#include "../src/constants.h"

void setUp(){}
void tearDown(){}

void assert_equal_complex_double_array(const double complex * expected,
                                        const double complex * actual,
                                        int length)
{
    int i;
    double complex elem1, elem2;
    char template[] = "Element %d %s Part";
    char message[40];
    for(i = 0; i < length; i++)
    {
        elem1 = *(actual + i);
        elem2 = *(expected + i);
        snprintf(message, 40, template, i, "Real");
        TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(crealf(elem2), crealf(elem1), message);
        snprintf(message, 40, template, i, "Imag");
        TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimagf(elem2), cimagf(elem1), message);

    }
}

void assert_equal_complex_float_array(const float complex * expected,
                                        const float complex * actual,
                                        int length)
{
    int i;
    float complex elem1, elem2;
    char template[] = "Element %d %s Part";
    char message[40];
    for(i = 0; i < length; i++)
    {
        elem1 = *(actual + i);
        elem2 = *(expected + i);
        snprintf(message, 40, template, i, "Real");
        TEST_ASSERT_EQUAL_FLOAT_MESSAGE(crealf(elem2), crealf(elem1), message);
        snprintf(message, 40, template, i, "Imag");
        TEST_ASSERT_EQUAL_FLOAT_MESSAGE(cimagf(elem2), cimagf(elem1), message);

    }
}


void test_io_read_until_char()
{
    char string[] = "csnaid8j3rnaw8caion4ojo\n7dh3ndkuay7ih^T&^bjuafhiji8hn*jk8h8\njbsadLIUNWIUkjds(Afundsfhxc\nashni87sdni3nnl38jkdil&b7y(Bkndsfkjljfdnmmlkmd\n";
    char filename[] = "data/test_file.txt";
    FILE * fp = fopen(filename, "w");
    fprintf(fp, string);
    fclose(fp);

    fp = fopen(filename, "r");
    io_read_until_char('(', fp);
    fseek(fp, 1, SEEK_CUR);
    char read_value = fgetc(fp);
    TEST_ASSERT_EQUAL_CHAR('A', read_value);

    io_read_until_char('(', fp);
    fseek(fp, 1, SEEK_CUR);
    read_value = fgetc(fp);
    TEST_ASSERT_EQUAL_CHAR('B', read_value);

    fclose(fp);
}

void test_io_read_array_complex()
{
    // int maxchar = MAXCHAR;
    // printf("MAXCHAR: %d\n", maxchar);
    CDTYPE expected[12] = { 0,     1.1,       2.2,
                            0.1*I, 1.1-0.1*I, 2.2-0.1*I,
                            0.2*I, 1.1+0.2*I, 2.2+0.2*I,
                            0.3*I, 1.1+0.3*I, 2.2+0.3*I};
    CDTYPE actual[12];
    io_read_array('C', 'C', actual, 3, 4, "data/test_matrix.dat");
    // utils_print_matrix(actual, 3, 4, 'C', 'C');
    assert_equal_complex_float_array(expected, actual, 12);
}

void test_io_readwrite_array_complex()
{
    // int maxchar = MAXCHAR;
    // printf("MAXCHAR: %d\n", maxchar);
    CDTYPE expected[12] = { 0,     1.1,       2.2,
                            0.1*I, 1.1-0.1*I, 2.2-0.1*I,
                            0.2*I, 1.1+0.2*I, 2.2+0.2*I,
                            0.3*I, 1.1+0.3*I, 2.2+0.3*I};
    CDTYPE actual[12];
    io_write_array_bin('C', expected, 3, 4, "data/test_matrix_bin.dat");
    io_read_array_bin('C', actual, 3, 4, "data/test_matrix_bin.dat");
    // utils_print_matrix(actual, 3, 4, 'C', 'C');
    assert_equal_complex_float_array(expected, actual, 12);
}


int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_io_read_until_char);
    RUN_TEST(test_io_read_array_complex);
    RUN_TEST(test_io_readwrite_array_complex);
    return(UNITY_END());
}