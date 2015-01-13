#ifndef RDTSC_H
#define RDTSC_H

inline unsigned long long rdtsc(void) {
  unsigned long long ret;
  __asm__ volatile ("rdtsc" : "=A" (ret));
  return ret;
}

#endif /* RDTSC_H */
