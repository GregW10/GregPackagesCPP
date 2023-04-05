#include "simsup.hpp"
#include <sys/shm.h>
#include <sys/wait.h>

#define NUM_PROCS 128

#ifdef __linux__
#define FILE1 "/home/mario_garwh/coms/iden"
#define FILE2 "/home/mario_garwh/coms/mass"
#else
#define FILE1 "/Users/gregorhartlwatters/MyCommands/iden"
#define FILE2 "/Users/gregorhartlwatters/MyCommands/mass"
#endif

long double b_rad = 69.157838586870785307914033523957l;
long double b_starting_mass = 1215368547.210631793481297791l;

long double bulk_density = 500; // kg / m^3

gtd::vector3D<long double> pos{-414139744.3484770l, 277323640.2236369l, -1231468367.968793l};
gtd::vector3D<long double> vel{3491.263406628809l, -6314.208154956334l, 11582.30048080498l};

uint64_t num = 727;
long double bounding_rad = 1600.0l;
long double restc_f = 1.0l;
long double min_cor = 0.5l;
long double max_cor = 0.95l;
long double ev_dt = 0.0625l/8;

long double dist_tol = 0.0l;

long double n_exp = 3;
long double d_scale = 1200;

uint64_t probing_iters = 100'000;

void func(FILE *fp, long double *pf_ptr, long double *rad_ptr, pid_t pid) {
    if (!fp || !pf_ptr || !rad_ptr)
        return;
    time_t _t;
    fprintf(fp, "Starting comet generation for process with pid %d at: %s\n", pid, ctime(&(_t = time(nullptr))));
    fflush(fp);
    auto [sys, pf, crad] =
            gtd::system<long double, long double, long double, false, false, 3, 0, 0, false>::
            random_comet<false>(pos, vel, num, bounding_rad, b_starting_mass, b_rad, restc_f, n_exp, d_scale,
                               gtd::sys::leapfrog_kdk, min_cor, max_cor, ev_dt, probing_iters, dist_tol);
    fprintf(fp, "Comet generation ended at: %s\n", ctime(&(_t = time(nullptr))));
    long double b_mass = gtd::adjust_bd(sys, b_rad, pf, bulk_density);
    fprintf(fp, "Comet effective radius: %.30Lf m\nBulk density: %.30Lf kg/m^3\nPacking fraction: %.30Lf\n",
            crad, (b_mass*num)/SPHERE_VOLUME(crad), pf);
    fprintf(fp, "Starting mass of each body: %.30Lf, final: %.30Lf\n", b_starting_mass, b_mass);
    fclose(fp);
    *pf_ptr = pf;
    *rad_ptr = crad;
    sys.to_nsys(gtd::String{}.append_back(pid).append_back(".nsys").c_str());
}

int main(int argc, char **argv) {
    static_assert((NUM_PROCS & (NUM_PROCS - 1)) == 0 && NUM_PROCS > 0);
    int shmid_pf = shmget(ftok(FILE1, 405), NUM_PROCS*sizeof(long double), IPC_CREAT | 0666); // 110110110 - rw for ugo
    if (shmid_pf == -1) {
        std::cerr << "shmget() for packing fraction failed.\n";
        perror("Error");
        return 1;
    }
    long double *const pf_ptr = (long double *) shmat(shmid_pf, nullptr, 0);
    if (pf_ptr == (void *) -1) {
        std::cerr << "shmat() for packing fraction failed.\n";
        perror("Error");
        return 1;
    }
    int shmid_rad = shmget(ftok(FILE2, 404), NUM_PROCS*sizeof(long double), IPC_CREAT | 0666);
    if (shmid_rad == -1) {
        std::cerr << "shmget() for radius failed.\n";
        perror("Error");
        return 1;
    }
    long double *rad_ptr = (long double *) shmat(shmid_rad, nullptr, 0);
    if (rad_ptr == (void *) -1) {
        std::cerr << "shmat() for radius failed.\n";
        perror("Error");
        return 1;
    }
    long double *pf = pf_ptr;
    const unsigned char n = (unsigned char) log2(NUM_PROCS);
    pid_t *pids = new pid_t[n];
    pid_t *pidptr = pids;
    for (unsigned char i = 1; i <= n; ++i) {
        if ((*pidptr = fork()) == -1) {
            std::cerr << "fork() error in process with pid: " << getpid() << std::endl;
            perror("Error that occurred");
            return 1;
        }
        if (!*pidptr++) {
            pf += 1 << (n - i);
            rad_ptr += 1 << (n - i);
        }
    }
    delete [] pids;
    unsigned int diff = (unsigned int) (pf - pf_ptr);
    sleep(diff ? 2*diff + 10 : 0);
    pid_t pid = getpid();
    FILE *fp = fopen(gtd::String{}.append_back(pid).append_back(".txt").c_str(), "w");
    func(fp, pf, rad_ptr, pid);
    if (pf != pf_ptr)
        return 0;
    while (wait(nullptr) > 0);
    printf("Average packing fraction and effective radius between %d random comets: %.30Lf & %.30Lf m\n",
           NUM_PROCS, gtd::mean_avg(pf, NUM_PROCS), gtd::mean_avg(rad_ptr, NUM_PROCS));
    shmdt(pf_ptr);
    shmctl(shmid_pf, IPC_RMID, nullptr);
    shmdt(rad_ptr);
    shmctl(shmid_rad, IPC_RMID, nullptr);
    return 0;
}
