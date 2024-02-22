#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HASHSIZE 211

/* KeyValuePair.key is an integer key
KeyValuePair.value is an integer value
KeyValuePair.next is a pointer to the next KeyValuePair whose key has the same hash value */
struct KeyValuePair {
    int key;
    int value;
    struct KeyValuePair *next;
};

/* Dict.hash_table[i] is a pointer to a KeyValuePair whose hashed key is i */
typedef struct {
    unsigned hash_size;
    struct KeyValuePair** hash_table;
} Dict;

unsigned hash(int k, unsigned hash_size) {
    /* Simple hash function. Modular arithmetic. */
    return k % hash_size;
}

struct KeyValuePair* lookup(Dict dict, int key) {
    /* Performs hash table lookup of a dictionary
    
    Inputs
    ------
    dict, pointer to dictionary
    key, lookup key
    
    Returns
    -------
    pointer to the required key-value pair */

    struct KeyValuePair* item;
    for (item=dict.hash_table[hash(key, dict.hash_size)]; item!=NULL; item=item->next) {
        if (key==item->key) {
            return item; // found
        }
    }
    return NULL; // not found
}

void print_dictionary(Dict dict) {
    /* Prints all key-value pairs of a dictionary in the format of key:value.
    
    Inputs
    ------
    dict, dictionary */
    for (int i=0; i<dict.hash_size; ++i) {
        struct KeyValuePair* item;
        for (item=dict.hash_table[i]; item!=NULL; item=item->next) {
            printf("%d:%d\n", item->key, item->value);
        }
    }
}


Dict counter(int *array, int length, int hash_size) {
    /* Returns a dictionary which counts the occurences of unique integers in array
    
    Inputs
    ------
    array, 1d array of ints
    length, length of array */

    // Initialize dictionary with hash_table pointers to NULL
    Dict dict;
    dict.hash_size = hash_size;
    dict.hash_table = malloc(dict.hash_size * sizeof(struct KeyValuePair*));
    for (int h=0; h<dict.hash_size; ++h) {dict.hash_table[h] = NULL;}

    for (int i=0; i<length; ++i) {
        struct KeyValuePair* item = lookup(dict, array[i]);

        // Install new key-value pair to dictionary if array[i] not found
        if (item == NULL) {
            unsigned hash_value = hash(array[i], dict.hash_size);
            item = malloc(sizeof(*item));
            item->key = array[i];
            item->value = 1;
            item->next = dict.hash_table[hash_value];
            dict.hash_table[hash_value] = item;
        } else {
            item->value++; // increment count by 1 when array[i] is found
        }
    }
    return dict;
}

int main() {
    int a[8] = {100, 100, 300, 200, 300, 1, 212, 212};
    Dict counts = counter(a, 8, HASHSIZE);
    print_dictionary(counts);
    free(counts.hash_table);
    printf("freed\n");
    return 0;
}