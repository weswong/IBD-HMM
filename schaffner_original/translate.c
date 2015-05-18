int base2index(char base_char) {
  int base;
  switch(base_char) {
  case 'A':
    base = 1;
    break;
  case 'C':
    base = 2;
    break;
  case 'G':
    base = 3;
    break;
  case 'T':
    base = 4;
    break;
  default: 
    base = 0;
  }
  return base;
}

char index2base(int index) {
  char base;
  switch(index) {
  case 1:
    base = 'A';
    break;
  case 2:
    base = 'C';
    break;
  case 3:
    base = 'G';
    break;
  case 4:
    base = 'T';
    break;
  default:
    base = 'N';
  }
  return base;
}
