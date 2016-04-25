#define NUM_FEATURES 6
#define INST_LEN (NUM_FEATURES * 8)
#define N_ULONG_IN_RULE 51
#define RULE_CASES 18 // was 17 with original short-circuiting rule matching
#define RULE_LEN (RULE_CASES * NUM_FEATURES * 8) // was 1632 based on 17 RULE_CASES and 12 Features
#define ONE_FREQ 6
#define MAX_NUM_CROSS_BREED 20
#define PROB_REMAIN 0.2
#define TAG 17
#define WHICH_SELECT 4 // tournament [1] default, fitness proportional [2], diversity guided tournament select [3]
#define TOURNAMENT_SIZE 5
#define RELEVANCE_START 0.25
#define RELEVANCE_END 0.25
#define RELEVANCE_INCREMENT 0.25
#define NUM_CYCLES 1
#define RANDOM_SEED 0 //1461439101
#define WHEN_FULL_DATA_MULTIPLY 1
#define WHEN_FULL_DATA_THRESHOLD 10
#define WHEN_FULL_DATA_START 100
#define WHICH_MIGRATION_DEFAULT 'r'
#define WHEN_MIGRATE_DEFAULT 50001
#define NUM_MIGRANTS_PER_ISLAND_DEFAULT 0
#define NUM_ISLANDS_DEFAULT 1
#define NUM_NEIGHBORS_DEFAULT 0
#define POP_SIZE_DEFAULT 400
#define NUM_TEST_CASES_TO_USE_DEFAULT 1000
#define MAX_GENERATIONS_DEFAULT 100
#define WHICH_FITNESS_DEFAULT 'h'
#define WHICH_CLASSIFY_DEFAULT 'h'
