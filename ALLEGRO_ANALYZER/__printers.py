from ___constants_names import PRGM_VERSION, PRGM_COOL_NAME, PRGM_END_CARD
import datetime

def print_start_msg():
    print(PRGM_COOL_NAME)
    print('Version:', PRGM_VERSION)
    start_time = datetime.datetime.now()
    print('Time program started:', start_time, '\n')
    return start_time

def print_end_msg(start_time):
    print(PRGM_END_CARD)
    end_time = datetime.datetime.now()
    print('Time program ended:', end_time, '\n')
    time_elapsed = end_time - start_time
    time_elapsed = divmod(time_elapsed.total_seconds(), 60)
    print('Total elapsed time: %f minutes, %f seconds.'%(time_elapsed[0], time_elapsed[1]))