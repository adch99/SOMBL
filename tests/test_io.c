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

void test_io_read_bin_d2f_real()
{
    double expected[10] = {0.85689756,  0.10860182,
                        1.58300039, -0.39497725,
                        0.13802642, -1.46729247,
                        -0.99387897, -1.87410812,
                        -0.07429885,  1.46899403};

    char filename[32] = "data/test_matrix_bin_d2f.dat";
    FILE * fp = io_safely_open_binary('w', filename);
    fwrite(expected, sizeof(double), 10, fp);
    fclose(fp);

    DTYPE actual[10];
    DTYPE diff;
    io_read_array_bin_d2f('R', actual, 5, 2, filename);
    int i;
    for(i = 0; i < 10; i++)
    {
        // printf("%le\t%e\n", expected[i], actual[i]);
        diff = fabsf((DTYPE) expected[i] - actual[i]);
        TEST_ASSERT(diff < 1e-12);
    }
}

void test_io_read_bin_d2f_complex()
{
    complex double expected[10] = {
        1.40486024-0.44140975*I, -1.57381065+0.37067508*I,
        0.22793477+1.8792268*I , -1.38351069-1.84556993*I,
       -0.13283254-1.38679858*I,  0.42847776-0.89848716*I,
        1.64092895+0.67027061*I, -0.26917277-1.33239699*I,
       -1.72383996+0.64206111*I, -0.2081726 -1.80638975*I};

    char filename[32] = "data/test_matrix_bin_d2f.dat";
    FILE * fp = io_safely_open_binary('w', filename);
    fwrite(expected, sizeof(complex double), 10, fp);
    fclose(fp);

    CDTYPE actual[10];
    DTYPE diff;
    io_read_array_bin_d2f('C', actual, 5, 2, filename);
    int i;
    for(i = 0; i < 10; i++)
    {
        // printf("%le\t%e\n", expected[i], actual[i]);
        diff = cabsf((CDTYPE) expected[i] - actual[i]);
        TEST_ASSERT(diff < 1e-12);
    }

}

void test_io_get_initial_vector()
{
    char filename[] = "data/testing_sample_initial_cond_L4.dat";
    DTYPE expected[32] = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,
                    0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1};
    DTYPE actual[32];
    int i;
    for(i = 0; i < 32; i++)
        actual[i] = 0;
    io_get_initial_cond_vector(actual, 'R', filename);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-12, expected, actual, 32);

    CDTYPE expectedc[32] = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,
                    0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1,
                    0, 0, 0, 0, 1, 1, 1};
    CDTYPE actualc[32];
    for(i = 0; i < 32; i++)
        actualc[i] = 0;
    io_get_initial_cond_vector(actualc, 'C', filename);
    assert_equal_complex_float_array(expectedc, actualc, 32);
}


int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_io_read_until_char);
    RUN_TEST(test_io_read_array_complex);
    RUN_TEST(test_io_readwrite_array_complex);
    RUN_TEST(test_io_read_bin_d2f_real);
    RUN_TEST(test_io_read_bin_d2f_complex);
    RUN_TEST(test_io_get_initial_vector);
    return(UNITY_END());
}