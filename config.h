#define NUM_FEATURES 6
#define INST_LEN (NUM_FEATURES * 8)
#define N_ULONG_IN_RULE 51
#define RULE_CASES 18 // was 17 with original short-circuiting rule matching
#define RULE_LEN (RULE_CASES * NUM_FEATURES * 8) // was 1632 based on 17 RULE_CASES and 12 Features
#define ONE_FREQ 6
#define MAX_NUM_CROSS_BREED 20
#define PROB_REMAIN 0.2
#define TAG 17
#define WHICH_SELECT 0 // tournament [1] default, fitness proportional [2], diversity guided tournament select [3]
#define TOURNAMENT_SIZE 5
#define RELAVANCE_WEIGHT 1.0
