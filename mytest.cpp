#include "vacdb.h"
#include <math.h>
#include <random>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

enum RANDOM {
    UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE
};

class Random {
public:
    Random(int min, int max, RANDOM type = UNIFORMINT, int mean = 50, int stdev = 20) : m_min(min), m_max(max),
                                                                                        m_type(type) {
        if (type == NORMAL) {
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor
            m_normdist = std::normal_distribution<>(mean, stdev);
        } else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min, max);
        } else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double) min, (double) max);
        } else { //the case of SHUFFLE to generate every number only once
            m_generator = std::mt19937(m_device());
        }
    }

    void setSeed(int seedNum) {
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    void getShuffle(vector<int> &array) {
        // the user program creates the vector param and passes here
        // here we populate the vector using m_min and m_max
        for (int i = m_min; i <= m_max; i++) {
            array.push_back(i);
        }
        shuffle(array.begin(), array.end(), m_generator);
    }

    void getShuffle(int array[]) {
        // the param array must be of the size (m_max-m_min+1)
        // the user program creates the array and pass it here
        vector<int> temp;
        for (int i = m_min; i <= m_max; i++) {
            temp.push_back(i);
        }
        std::shuffle(temp.begin(), temp.end(), m_generator);
        vector<int>::iterator it;
        int i = 0;
        for (it = temp.begin(); it != temp.end(); it++) {
            array[i] = *it;
            i++;
        }
    }

    int getRandNum() {
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if (m_type == NORMAL) {
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while (result < m_min || result > m_max)
                result = m_normdist(m_generator);
        } else if (m_type == UNIFORMINT) {
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum() {
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result * 100.0) / 100.0;
        return result;
    }

    string getRandString(int size) {
        // the parameter size specifies the length of string we ask for
        // to use ASCII char the number range in constructor must be set to 97 - 122
        // and the Random type must be UNIFORMINT (it is default in constructor)
        string output = "";
        for (int i = 0; i < size; i++) {
            output = output + (char) getRandNum();
        }
        return output;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution
};


int hashFunction(string str);


class Tester {
public:
    bool testInsertionNormal();

    bool testFindError();
    bool testFindWithNonCollidingKeys();
    bool testFindWithCollidingKeys();

    bool testRemoveWithNonCollidingKeys();
    bool testRemoveWithCollidingKeys();

    bool testRehashAfterInsertion();
    bool testRehashCompletionFromLoadFactor();
    bool testRehashAfterRemoval();
    bool testRehashCompletionFromDeletedRatio();

private:

};


bool Tester::testInsertionNormal() {
    //TODO: Test the insertion operation in the hash table. The following presents a sample algorithm to test the normal insertion operation:
    // There are some non-colliding data points in the hash table.
    // Insert multiple non-colliding keys.
    // Check whether they are inserted in the correct bucket (correct index).
    // Check whether the data size changes correctly.

    return false;
}

bool Tester::testFindError() {
    //TODO: Test the find operation (getPatient(...) function) for an error case, the Patient object does not exist in the database.

    return false;
}

bool Tester::testFindWithNonCollidingKeys() {
    //TODO: Test the find operation (getPatient(...) function) with several non-colliding keys.

    return false;
}

bool Tester::testFindWithCollidingKeys() {
    //TODO: Test the find operation (getPatient(...) function) with several colliding keys without triggering a rehash. This also tests whether the insertion works correctly with colliding data.

    return false;
}

bool Tester::testRemoveWithNonCollidingKeys() {
    //TODO: Test the remove operation with a few non-colliding keys.

    return false;
}

bool Tester::testRemoveWithCollidingKeys() {
    //TODO: Test the remove operation with a number of colliding keys without triggering a rehash.

    return false;
}

bool Tester::testRehashAfterInsertion() {
    //TODO: Test the rehashing is triggered after a descent number of data insertion.

    return false;
}

bool Tester::testRehashCompletionFromLoadFactor() {
    //TODO: Test the rehash completion after triggering rehash due to load factor, i.e. all live data is transferred to the new table and the old table is removed.

    return false;
}

bool Tester::testRehashAfterRemoval() {
    //TODO: Test the rehashing is triggered after a descent number of data removal.

    return false;
}

bool Tester::testRehashCompletionFromDeletedRatio() {
    //TODO: Test the rehash completion after triggering rehash due to delete ratio, i.e. all live data is transferred to the new table and the old table is removed.

    return false;
}


int main() {
    Tester tester;

    return 0;
}


int hashFunction(string id) {
    const int prime = 31;
    int result = 0;
    for (int i = 0; i < id.length(); i++) {
        result += id[i] * pow(prime, i);
    }
    return result;
}