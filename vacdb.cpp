// CMSC 341 - Spring 2024 - Project 4
#include "vacdb.h"
VacDB::VacDB(int size, hash_fn hash, prob_t probing = DEFPOLCY){
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
    m_currentTable = new Patient*[m_currentCap]();

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

VacDB::~VacDB(){
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

void VacDB::changeProbPolicy(prob_t policy){
    m_newPolicy = policy;
}

bool VacDB::insert(Patient patient){
    //check whether a duplicate patient already exists
    for (int i = 0; i < m_currentCap; i++) {
        if (*m_currentTable[i] == patient) {
            return false;
        }
    }

    //calculate index for insertion
    unsigned int index = m_hash(patient.getKey()) % m_currentCap;

    //hash collisions are resolved using the probing policy
    index = probe(index, patient.getKey());

    //if a patient already exists at calculated index, deallocate it
    delete m_currentTable[index];

    //insert patient at calculated index if serial number is valid
    if (patient.getSerial() >= MINID && patient.getSerial() <= MAXID) {
        //allocate memory for Patient object
        Patient* newPatient = new Patient(patient.getKey(), patient.getSerial(), patient.getUsed());

        //insert and update current number of entries
        m_currentTable[index] = newPatient;
        m_currentSize++;

        //if the load factor exceeds 50% after an insertion, rehash to a new hash table
        if (lambda() > 0.5) {
            rehash();
        }

        return true;
    }
    //if serial number is invalid, do not insert and return false
    return false;
}

unsigned int VacDB::probe(unsigned int index, string key) {
    int i = 0;
    while (m_currentTable[index] != nullptr) {
        switch (m_currProbing) {
            case LINEAR:
                index = (index + 1) % m_currentCap;
                break;
            case QUADRATIC:
                index = (index + i * i) % m_currentCap;
                break;
            case DOUBLEHASH:
                index = ((m_hash(key) % m_currentCap) + i * (11 - (m_hash(key) % 11))) % m_currentCap;
                break;
        }
        i++;
    }
    return index;
}

void VacDB::rehash() {
    //if on initial rehash
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

        m_currentTable = new Patient*[m_currentCap]();

        //calculate how many data points will be transferred in each rehash
        int percentToTransfer = floor(0.25 * numDataPoints);
    }

    //
    int numTransferred = 0;

    //
    for (; m_transferIndex < m_oldCap || numTransferred >= percentToTransfer; m_transferIndex++) {
        //only transfer live data to new table
        if (!m_oldTable[m_transferIndex]->getUsed()) {
            //calculate index for insertion
            unsigned int newIndex = m_hash(m_oldTable[m_transferIndex]->getKey()) % m_currentCap;

            //hash collisions are resolved using the probing policy
            newIndex = probe(newIndex, m_oldTable[m_transferIndex]->getKey());

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

bool VacDB::remove(Patient patient){
    //search for patient in the current table
    for (int i = 0; i < m_currentCap; i++) {
        if (*m_currentTable[i] == patient) {
            //if found, mark patient as deleted (soft-delete) and increment deleted count
            m_currentTable[i]->setUsed(false);
            m_currNumDeleted++;

            //if the deleted ratio exceeds 80% after a deletion, rehash
            if (deletedRatio() > 0.8) {
                rehash();
            }

            return true;
        }
    }

    //search for patient in the old table
    for (int i = 0; i < m_oldCap; i++) {
        if (*m_oldTable[i] == patient) {
            //if found, mark patient as deleted (soft-delete) and increment deleted count
            m_oldTable[i]->setUsed(false);
            m_oldNumDeleted++;

            //if the deleted ratio exceeds 80% after a deletion, rehash
            if (deletedRatio() > 0.8) {
                rehash();
            }

            return true;
        }
    }

    //if patient is not found, return false
    return false;
}

const Patient VacDB::getPatient(string name, int serial) const{
    //return the Patient object with the passed-in name and the vaccine serial number in the database
    for (int i = 0; i < m_currentCap; i++) {
        if (m_currentTable[i]->getKey() == name && m_currentTable[i]->getSerial() == serial) {
            return *m_currentTable[i];
        }
    }
    for (int i = 0; i < m_oldCap; i++) {
        if (m_oldTable[i]->getKey() == name && m_oldTable[i]->getSerial() == serial) {
            return *m_oldTable[i];
        }
    }

    //if object is not found, return empty object
    return Patient();
}

bool VacDB::updateSerialNumber(Patient patient, int serial){
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

bool VacDB::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VacDB::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) {
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0)
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

ostream& operator<<(ostream& sout, const Patient* patient ) {
    if ((patient != nullptr) && !(patient->getKey().empty()))
        sout << patient->getKey() << " (" << patient->getSerial() << ", "<< patient->getUsed() <<  ")";
    else
        sout << "";
    return sout;
}

bool operator==(const Patient& lhs, const Patient& rhs){
    // since the uniqueness of an object is defined by name and serial number
    // the equality operator considers only those two criteria
    return ((lhs.getKey() == rhs.getKey()) && (lhs.getSerial() == rhs.getSerial()));
}

bool Patient::operator==(const Patient* & rhs){
    // since the uniqueness of an object is defined by name and serial number
    // the equality operator considers only those two criteria
    return ((getKey() == rhs->getKey()) && (getSerial() == rhs->getSerial()));
}