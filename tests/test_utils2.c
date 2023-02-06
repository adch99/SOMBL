#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../extern/unity/unity.h"
#include "../src/utils/utils.h"
#include "../src/io/io.h"

void setUp() {}
void tearDown() {}

void assert_equal_complex_float_array(float complex * expected,
                                        float complex * actual,
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
        snprintf(message, 40, template, i, "Complex");
        TEST_ASSERT_EQUAL_FLOAT_MESSAGE(cimagf(elem2), cimagf(elem1), message);

    }
}

void test_utils_construct_data_vs_dist()
{
    
}

void test_utils_add_to_matrix_real()
{
    float matrix1[8] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1};
    float matrix2[8] = {10.2, 20.2, 30.2, 40.2, 50.2, 60.2, 70.2, 80.2};
    float expected[8] = {11.3, 22.3, 33.3, 44.3, 55.3, 66.3, 77.3, 88.3};

    utils_add_to_matrix_real(matrix1, matrix2, 4, 2);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected, matrix1, 8);
}

void test_utils_add_to_matrix_complex()
{
    float complex matrix1[8] = {1.1 + I*0.2, 2.1 + I*0.2, 3.1 + I*0.2, 4.1 + I*0.2,
                                5.1 + I*0.2, 6.1 + I*0.2, 7.1 + I*0.2, 8.1 + I*0.2};
    float complex matrix2[8] = {10.2 + I*1.1, 20.2 + I*2.1, 30.2 + I*3.1, 40.2 + I*4.1,
                                50.2 + I*5.1, 60.2 + I*6.1, 70.2 + I*7.1, 80.2 + I*8.1};
    float complex expected[8] = {11.3 + I*1.3, 22.3 + I*2.3, 33.3 + I*3.3, 44.3 + I*4.3,
                                55.3 + I*5.3, 66.3 + I*6.3, 77.3 + I*7.3, 88.3 + I*8.3};

    utils_add_to_matrix_complex(matrix1, matrix2, 4, 2);
    assert_equal_complex_float_array(expected, matrix1, 8);
    // int i, j;
    // CDTYPE elem1; elem2;

}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_utils_add_to_matrix_real);
    RUN_TEST(test_utils_add_to_matrix_complex);
    return(UNITY_END());
}