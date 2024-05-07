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


unsigned int hashFunction(string str);

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
    vector<Patient> insertMultiplePatients(VacDB &vaccineDatabase);

    bool probe(unsigned int &index, string key, int serial, bool isCurrentTable, VacDB &vaccineDatabase) const;
};


vector<Patient> Tester::insertMultiplePatients(VacDB &vaccineDatabase) {
    Random randKeyObject(97, 122);
    Random randSerialObject(MINID, MAXID);

    int patientSize = 25;

    vector<Patient> patientVector;

    for (int i = 0; i < patientSize; i++) {
        string randKey = randKeyObject.getRandString(10);
        int randSerial = randSerialObject.getRandNum();

        Patient patient(randKey, randSerial, true);
        vaccineDatabase.insert(patient);

        patientVector.push_back(patient);
    }
    return patientVector;
}

bool Tester::probe(unsigned int &index, string key, int serial, bool isCurrentTable, VacDB &vaccineDatabase) const {
    //declare required variables
    Patient **hashTable;
    int capacity;
    prob_t probingPolicy;

    bool softDeleteFound = false;
    unsigned int firstSoftDeletedIndex = 0;

    //table in question changes based on boolean passed in
    if (isCurrentTable) {
        hashTable = vaccineDatabase.m_currentTable;
        capacity = vaccineDatabase.m_currentCap;
        probingPolicy = vaccineDatabase.m_currProbing;
    } else {
        hashTable = vaccineDatabase.m_oldTable;
        capacity = vaccineDatabase.m_oldCap;
        probingPolicy = vaccineDatabase.m_oldProbing;
    }

    //if table is empty
    if (hashTable == nullptr) {
        return false;
    }

    //calculate initial index
    index = vaccineDatabase.m_hash(key) % capacity;

    //loop through table until an empty slot or match is found
    for (int i = 1; hashTable[index] != nullptr; i++) {
        //save first soft-deleted index in new table
        if (hashTable == vaccineDatabase.m_currentTable && !hashTable[index]->getUsed() && !softDeleteFound) {
            softDeleteFound = true;
            firstSoftDeletedIndex = index;
        }
        //if match found, return true
        if (hashTable[index]->getKey() == key && hashTable[index]->getSerial() == serial && hashTable[index]->m_used) {
            return true;
        }
        //increment index via probing policy
        switch (probingPolicy) {
            case LINEAR:
                index = (index + 1) % capacity;
                break;
            case QUADRATIC:
                index = (index + i * i) % capacity;
                break;
            case DOUBLEHASH:
                index = ((vaccineDatabase.m_hash(key) % capacity) + i * (11 - (vaccineDatabase.m_hash(key) % 11))) %
                        capacity;
                break;
        }
    }

    //if match not found but a soft-delete is found, index should change
    //soft-deleted index has priority over empty index
    if (softDeleteFound) {
        index = firstSoftDeletedIndex;
    }

    return false;
}

bool Tester::testInsertionNormal() {
    //insert multiple non-colliding data points into the hash table
    VacDB vaccineDatabase(MINPRIME, hashFunction, DOUBLEHASH);
    vector<Patient> patientVector = insertMultiplePatients(vaccineDatabase);

    //check whether they are inserted in the correct bucket (correct index)
    for (unsigned int i = 0; i < patientVector.capacity(); i++) {
        if (!probe(i, patientVector[i].getKey(), patientVector[i].getSerial(), true, vaccineDatabase)) {
            return false;
        }
    }

    //check whether the data size changes correctly
    if (vaccineDatabase.m_currentSize == 25) {
        return true;
    }
    return false;
}

bool Tester::testFindError() {
    VacDB vaccineDatabase(MINPRIME, hashFunction, DOUBLEHASH);

    //insert one patient
    Patient insertedPatient("Ymir", 1999, true);
    vaccineDatabase.insert(insertedPatient);

    //search for a patient that does not exist in the database
    Patient foundPatient = vaccineDatabase.getPatient("Maria", 2012);

    //if an empty object is returned, test passes
    if (foundPatient.getKey().empty() && foundPatient.getSerial() == 0 && !foundPatient.getUsed()) {
        return true;
    }
    return false;
}

bool Tester::testFindWithNonCollidingKeys() {
    //insert multiple non-colliding data points into the hash table
    VacDB vaccineDatabase(MINPRIME, hashFunction, DOUBLEHASH);
    vector<Patient> patientVector = insertMultiplePatients(vaccineDatabase);

    //search for an inserted patient
    Patient foundPatient = vaccineDatabase.getPatient(patientVector[0].getKey(), patientVector[0].getSerial());

    //if the found patient's data matches that of the one we were looking for, test passes
    if (foundPatient.getKey() == patientVector[0].getKey() &&
        foundPatient.getSerial() == patientVector[0].getSerial()) {
        return true;
    }
    return false;
}

bool Tester::testFindWithCollidingKeys() {
    VacDB vaccineDatabase(MINPRIME, hashFunction, DOUBLEHASH);

    Random randSerialObject(MINID, MAXID);

    int patientSize = 25;

    vector<Patient> patientVector;

    //insert colliding patients
    for (int i = 0; i < patientSize; i++) {
        int randSerial = randSerialObject.getRandNum();

        Patient patient("Ymir", randSerial, true);
        vaccineDatabase.insert(patient);

        patientVector.push_back(patient);
    }

    //search for an inserted patient
    Patient foundPatient = vaccineDatabase.getPatient("Ymir", patientVector[0].getSerial());

    //if the found patient's data matches that of the one we were looking for, test passes
    if (foundPatient.getKey() == "Ymir" && foundPatient.getSerial() == patientVector[0].getSerial()) {
        return true;
    }
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

    cout << "\nTesting insert (normal) - check insertion with multiple non-colliding keys:" << endl;
    if (tester.testInsertionNormal()) {
        cout << "\tTest passed!" << endl;
    } else {
        cout << "\t***Test failed!***" << endl;
    }

    cout << "\nTesting getPatient (error) - ___:" << endl;
    if (tester.testFindError()) {
        cout << "\tTest passed!" << endl;
    } else {
        cout << "\t***Test failed!***" << endl;
    }
    cout << "Testing getPatient (non-colliding data) - ___:" << endl;
    if (tester.testFindWithNonCollidingKeys()) {
        cout << "\tTest passed!" << endl;
    } else {
        cout << "\t***Test failed!***" << endl;
    }
    cout << "Testing getPatient (colliding data) - ___:" << endl;
    if (tester.testFindWithCollidingKeys()) {
        cout << "\tTest passed!" << endl;
    } else {
        cout << "\t***Test failed!***" << endl;
    }

    return 0;
}


unsigned int hashFunction(string id) {
    const int prime = 31;
    int result = 0;
    for (unsigned int i = 0; i < id.length(); i++) {
        result += id[i] * pow(prime, i);
    }
    return result;
}