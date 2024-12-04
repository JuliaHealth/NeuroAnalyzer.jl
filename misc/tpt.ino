// Two-point Pinch Test (TPT)

#include <Wire.h>
#include "MMA7660.h"

MMA7660 accelemeter;

void setup()
{
	accelemeter.init();  
	Serial.begin(19200);
}

void loop()
{
  // position
	int8_t x, y, z;
  // acceleration
	float ax, ay, az;
  float g = 9.8066;
  // timing
  unsigned long t1, t2, t, d;
  // sampling frequency: 50 Hz
  int fs = 50;

  // microseconds
  t1 = micros();

  accelemeter.getXYZ(&x, &y, &z);
	accelemeter.getAcceleration(&ax, &ay, &az);

  // convert from g to m/s2
  ax *= g;
  ay *= g;
  az *= g;

  Serial.print("tpt: ");
  Serial.print(x); 
  Serial.print(" ");
  Serial.print(y);
  Serial.print(" ");
  Serial.print(z);
  Serial.print(" ");
	Serial.print(ax, 3);
  Serial.print(" ");
	Serial.print(ay, 3);
  Serial.print(" ");
	Serial.println(az, 3);

  // microseconds
  t2 = micros();

  t = (t2 - t1) / 1000.0;
  d = (1000.0 / fs) - t;
	delay(d);
}


