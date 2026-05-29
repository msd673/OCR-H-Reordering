#include "eval.h"

int perf_event_open_wrapper(struct perf_event_attr* hw_event, pid_t pid,
                            int cpu, int group_fd, unsigned long flags, const string& event_name) {
    int fd = syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
    if (fd == -1) {
        cerr << "perf_event_open (" << event_name << ") failed: " << strerror(errno) 
             << " (Try: sudo sh -c \"echo 0 > /proc/sys/kernel/perf_event_paranoid\" or check permissions)"
             << endl;
    }
    return fd;
}

CachePerfResult measure_cache_perf(function<void()> workload) {
    CachePerfResult result = {};
    
    int fd_l1_ref = -1;
    int fd_l1_miss = -1;
    int fd_llc_access = -1;
    int fd_llc_miss = -1;

    struct perf_event_attr pe_l1;
    memset(&pe_l1, 0, sizeof(pe_l1));
    pe_l1.size = sizeof(pe_l1);
    pe_l1.type = PERF_TYPE_HW_CACHE;
    pe_l1.config = (PERF_COUNT_HW_CACHE_L1D << 0) |
                   (PERF_COUNT_HW_CACHE_OP_READ << 8) |
                   (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16);
    pe_l1.disabled = 1;
    pe_l1.exclude_kernel = 1;
    pe_l1.exclude_hv = 1;
    pe_l1.read_format = PERF_FORMAT_GROUP;

    fd_l1_ref = perf_event_open_wrapper(&pe_l1, 0, -1, -1, 0, "L1_REF");
    if (fd_l1_ref == -1) {
        return result;
    }

    pe_l1.config = (PERF_COUNT_HW_CACHE_L1D << 0) |
                   (PERF_COUNT_HW_CACHE_OP_READ << 8) |
                   (PERF_COUNT_HW_CACHE_RESULT_MISS << 16);
    
    fd_l1_miss = perf_event_open_wrapper(&pe_l1, 0, -1, fd_l1_ref, 0, "L1_MISS");
    if (fd_l1_miss == -1) {
        if (fd_l1_ref != -1) close(fd_l1_ref);
        return result;
    }

    struct perf_event_attr pe_llc;
    memset(&pe_llc, 0, sizeof(pe_llc));
    pe_llc.size = sizeof(pe_llc);
    pe_llc.type = PERF_TYPE_HW_CACHE;
    pe_llc.config = (PERF_COUNT_HW_CACHE_LL << 0) |
                     (PERF_COUNT_HW_CACHE_OP_READ << 8) |
                     (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16);
    pe_llc.disabled = 1;
    pe_llc.read_format = PERF_FORMAT_GROUP;
    
    fd_llc_access = perf_event_open_wrapper(&pe_llc, 0, -1, -1, 0, "LLC_ACCESS");
    if (fd_llc_access == -1) {
        if (fd_l1_ref != -1) close(fd_l1_ref);
        if (fd_l1_miss != -1) close(fd_l1_miss);
        return result;
    }

    pe_llc.config = (PERF_COUNT_HW_CACHE_LL << 0) |
                     (PERF_COUNT_HW_CACHE_OP_READ << 8) |
                     (PERF_COUNT_HW_CACHE_RESULT_MISS << 16);
    
    fd_llc_miss = perf_event_open_wrapper(&pe_llc, 0, -1, fd_llc_access, 0, "LLC_MISS");
    if (fd_llc_miss == -1) {
        if (fd_l1_ref != -1) close(fd_l1_ref);
        if (fd_l1_miss != -1) close(fd_l1_miss);
        if (fd_llc_access != -1) close(fd_llc_access);
        return result;
    }

    if (ioctl(fd_l1_ref, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == -1) {
        cerr << "Enable L1 group failed: " << strerror(errno) << endl;
    }
    
    if (ioctl(fd_llc_access, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == -1) {
        cerr << "Enable LLC group failed: " << strerror(errno) << endl;
    }

    workload();

    long long values[3];
    ssize_t bytes_read_l1 = read(fd_l1_ref, values, sizeof(values));
    if (bytes_read_l1 == sizeof(values)) {
        if (values[0] >= 2) { 
            result.l1_ref = values[1];
            result.l1_miss = values[2];
        } else {
            cerr << "Error: L1 group read, but expected at least 2 events in group, got " << values[0] << endl;
        }
    } else {
        cerr << "Error: Read L1 group failed. Bytes read: " << bytes_read_l1 << ", expected: " << sizeof(values) << ". System error: " << strerror(errno) << endl;
    }

    memset(values, 0, sizeof(values));
    ssize_t bytes_read_llc = read(fd_llc_access, values, sizeof(values));
    if (bytes_read_llc == sizeof(values)) {
        if (values[0] >= 2) { 
            result.llc_ref = values[1];
            result.llc_miss = values[2];
            result.cache_miss_total = values[2]; 
        } else {
            cerr << "Error: LLC group read, but expected at least 2 events in group, got " << values[0] << endl;
        }
    } else {
        cerr << "Error: Read LLC group failed. Bytes read: " << bytes_read_llc << ", expected: " << sizeof(values) << ". System error: " << strerror(errno) << endl;
    }

    if (fd_l1_ref != -1) ioctl(fd_l1_ref, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
    if (fd_llc_access != -1) ioctl(fd_llc_access, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);

    if (fd_l1_ref != -1) close(fd_l1_ref);
    if (fd_l1_miss != -1) close(fd_l1_miss);
    if (fd_llc_access != -1) close(fd_llc_access);
    if (fd_llc_miss != -1) close(fd_llc_miss);
    
    return result;
}

double measure_avg_time(int repeat, function<void()> func) {
    using namespace std::chrono;
    long long total_time = 0;

    for (int i = 0; i < repeat; ++i) {
        auto t1 = high_resolution_clock::now();
        func();
        auto t2 = high_resolution_clock::now();
        total_time += duration_cast<microseconds>(t2 - t1).count();
    }
    return total_time / (repeat * 1000.0);
}