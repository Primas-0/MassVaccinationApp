// CMSC 341 - Spring 2024 - Project 4
#include "vacdb.h"

VacDB::VacDB(int size, hash_fn hash, prob_t probing = DEFPOLCY) {
    //set table size after validation
    if (size < MINPRIME) {
        m_currentCap = MINPRIME;
    } else if (size > MAXPRIME) {
        m_currentCap = MAXPRIME;
    } else if (!isPrime(size)) {
        m_currentCap = findNextPrime(size);
    } else {
        m_currentCap = size;
    }

    //initialize current member variables
    m_hash = hash;
    m_newPolicy = probing;

    //create memory for the current table
    m_currentTable = new Patient *[m_currentCap]();

    //initialize all other member variables
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currProbing = probing;
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_oldProbing = probing;
    m_transferIndex = -1;
}

VacDB::~VacDB() {
    //deallocate current table and pointers
    for (int i = 0; i < m_currentCap; i++) {
        delete m_currentTable[i];
    }
    delete[] m_currentTable;

    //deallocate old table and pointers
    for (int i = 0; i < m_oldCap; i++) {
        delete m_oldTable[i];
    }
    delete[] m_oldTable;
}

void VacDB::changeProbPolicy(prob_t policy) {
    m_newPolicy = policy;
}

bool VacDB::insert(Patient patient) {
    //calculate index for insertion
    unsigned int index = m_hash(patient.getKey()) % m_currentCap;

    bool insertSuccessFlag = false;

    //hash collisions are resolved using the probing policy
    probe(index, patient.getKey(), patient.getSerial(), true);

    //if a patient already exists at calculated index, deallocate it
    delete m_currentTable[index];

    //insert patient at calculated index if there is no duplicate and serial number is valid
    if (!probe(index, patient.getKey(), patient.getSerial(), true) &&
        patient.getSerial() >= MINID && patient.getSerial() <= MAXID) {
        //allocate memory for Patient object
        Patient *newPatient = new Patient(patient.getKey(), patient.getSerial(), patient.getUsed());

        //insert and update current number of entries
        m_currentTable[index] = newPatient;
        m_currentSize++;

        //flag becomes true after insertion
        insertSuccessFlag = true;
    }

    //regardless of the output of the insert,
    //if the load factor exceeds 50% after an insertion, rehash (or if rehash is already in progress, continue)
    if (lambda() > 0.5 || m_transferIndex != -1) {
        rehash();
    }

    //if serial number is invalid, do not insert and return false
    return insertSuccessFlag;
}

bool VacDB::probe(unsigned int &index, string key, int serial, bool isCurrentTable) const {
    Patient **hashTable;
    int capacity;
    prob_t probingPolicy;

    bool softDeleteFound = false;
    unsigned int firstSoftDeletedIndex = 0;

    if (isCurrentTable) {
        hashTable = m_currentTable;
        capacity = m_currentCap;
        probingPolicy = m_currProbing;
    } else {
        hashTable = m_oldTable;
        capacity = m_oldCap;
        probingPolicy = m_oldProbing;
    }

    //loop through table until an empty slot or match is found
    for (int i = 0; hashTable[index] != nullptr; i++) {
        //save first soft-deleted index in new table
        if (!m_currentTable[index]->getUsed() && !softDeleteFound) {
            softDeleteFound = true;
            firstSoftDeletedIndex = index;
        }
        //if match found, return true
        if (hashTable[index]->getKey() == key && hashTable[index]->getSerial() == serial) {
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
                index = ((m_hash(key) % capacity) + i * (11 - (m_hash(key) % 11))) % capacity;
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

void VacDB::rehash() {
    //initial rehash setup
    if (m_transferIndex == -1) {
        //move current table to old table
        m_oldTable = m_currentTable;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_oldNumDeleted = m_currNumDeleted;
        m_oldProbing = m_currProbing;

        //clear current table member variables to use as the "new" table
        m_currentTable = nullptr;
        m_currentCap = 0;
        m_currentSize = 0;
        m_currNumDeleted = 0;
        m_currProbing = m_newPolicy;

        //calculate new table capacity
        int numDataPoints = m_oldSize - m_oldNumDeleted;
        m_currentCap = findNextPrime(4 * numDataPoints);

        m_currentTable = new Patient *[m_currentCap]();
    }

    //calculate how many data points will be transferred in each rehash
    int percentToTransfer = floor(0.25 * m_oldCap);

    int numTransferred = 0;

    //traverse through old table until it reaches the end and scan limit reached
    for (; m_transferIndex < m_oldCap && numTransferred < percentToTransfer; m_transferIndex++) {
        //only transfer live data to new table
        if (!m_oldTable[m_transferIndex]->getUsed()) {
            //calculate index for insertion
            unsigned int newIndex = m_hash(m_oldTable[m_transferIndex]->getKey()) % m_currentCap;

            //hash collisions are resolved using the probing policy
            probe(newIndex, m_oldTable[m_transferIndex]->getKey(), m_oldTable[m_transferIndex]->getSerial(), false);

            delete m_currentTable[newIndex];

            //allocate memory for Patient object
            Patient *newPatient = new Patient(m_oldTable[m_transferIndex]->getKey(),
                                              m_oldTable[m_transferIndex]->getSerial(),
                                              m_oldTable[m_transferIndex]->getUsed());

            //insert and update current number of entries
            m_currentTable[newIndex] = newPatient;
            m_currentSize++;

            numTransferred++;

            //soft-delete data from old table after inserting into new table
            m_oldTable[m_transferIndex]->setUsed(false);
        }
    }

    //remove old table and deallocate its memory once all rehashing is complete
    if (m_transferIndex == m_oldCap) {
        for (int i = 0; i < m_oldCap; i++) {
            delete m_oldTable[i];
        }
        delete[] m_oldTable;
    }
}

bool VacDB::remove(Patient patient) {
    bool removeSuccessFlag = false;

    unsigned int index = 0;

    //if found in either table, mark patient as deleted (soft-delete) and set success flag to true
    if (probe(index, patient.getKey(), patient.getSerial(), true)) {
        m_currentTable[index]->setUsed(false);
        m_currNumDeleted++;
        removeSuccessFlag = true;
    } else if (probe(index, patient.getKey(), patient.getSerial(), false)) {
        m_oldTable[index]->setUsed(false);
        m_oldNumDeleted++;
        removeSuccessFlag = true;
    }

    //regardless of the output of the remove,
    //if the deleted ratio exceeds 80% after a deletion, rehash (or if rehash is already in progress, continue)
    if (deletedRatio() > 0.8 || m_transferIndex != -1) {
        rehash();
    }

    //if patient is not found, return false
    return removeSuccessFlag;
}

const Patient VacDB::getPatient(string name, int serial) const {
    unsigned int index = 0;

    //return the Patient object with the passed-in name and the vaccine serial number in the database
    if (probe(index, name, serial, true)) {
        return *m_currentTable[index];
    } else if (probe(index, name, serial, false)) {
        return *m_oldTable[index];
    }

    //if object is not found, return empty object
    return Patient();
}

bool VacDB::updateSerialNumber(Patient patient, int serial) {
    //call getPatient to search the database
    Patient foundPatient = getPatient(patient.getKey(), patient.getSerial());

    //if the returned object is empty, return false
    if (foundPatient.getKey().empty() && foundPatient.getSerial() == 0 && !foundPatient.getUsed()) {
        return false;
    }

    //otherwise, update its serial number and return true
    foundPatient.setSerial(serial);
    return true;
}

float VacDB::lambda() const {
    //return load factor of current hash table
    return float(m_currentSize) / float(m_currentCap);
}

float VacDB::deletedRatio() const {
    //return the ratio of the deleted buckets to the total number of occupied buckets
    return float(m_currNumDeleted) / float(m_currentSize);
}

void VacDB::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool VacDB::isPrime(int number) {
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VacDB::findNextPrime(int current) {
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME - 1;
    for (int i = current; i < MAXPRIME; i++) {
        for (int j = 2; j * j <= i; j++) {
            if (i % j == 0)
                break;
            else if (j + 1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

ostream &operator<<(ostream &sout, const Patient *patient) {
    if ((patient != nullptr) && !(patient->getKey().empty()))
        sout << patient->getKey() << " (" << patient->getSerial() << ", " << patient->getUsed() << ")";
    else
        sout << "";
    return sout;
}

bool operator==(const Patient &lhs, const Patient &rhs) {
    // since the uniqueness of an object is defined by name and serial number
    // the equality operator considers only those two criteria
    return ((lhs.getKey() == rhs.getKey()) && (lhs.getSerial() == rhs.getSerial()));
}

bool Patient::operator==(const Patient *&rhs) {
    // since the uniqueness of an object is defined by name and serial number
    // the equality operator considers only those two criteria
    return ((getKey() == rhs->getKey()) && (getSerial() == rhs->getSerial()));
}