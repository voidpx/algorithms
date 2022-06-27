#ifndef CRC_H_
#define CRC_H_

#include <stddef.h>

#define POLY 0xEDB88320 /* 0x04C11DB7 bits reversed */

unsigned int calc_file_crc(char *file);
unsigned int calc_crc(char *data, size_t size);

#endif /* CRC_H_ */
