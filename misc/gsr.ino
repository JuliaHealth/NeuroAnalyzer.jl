const int gsr_sensor = A0;
double gsr_sensor_value = 0;
double gsr_average = 0;
double impedance = 0;

int calibration = 270;  // this may differ for different sensors, see https://wiki.seeedstudio.com/Grove-GSR_Sensor/ for details
int n = 10;             // length of moving average filter
long sum;
int fs = 20;            // 20 ms per loop, so sampling rate is 50 Hz
unsigned long t1;
unsigned long t2;

void setup(){
  Serial.begin(19200);
}

void loop(){
  t1 = millis();

  // record data and average the 10 measurements to remove the glitch
  sum = 0;
  for(int i = 0; i < n; i++) {
    gsr_sensor_value = analogRead(gsr_sensor);
    sum += gsr_sensor_value;
  }
  gsr_average = sum / n;

  // caclulate impedance
  impedance = 1000000 * (1 / (((1024 + 2 * gsr_average) * 10000) / (calibration - gsr_average))); // convert to uSiemens
  if (impedance < 0) {
    impedance = 0;
  }
  if (isnan(impedance)) {
    impedance = 0;
  }

  // send data
  Serial.print("gsr:");
  Serial.println(impedance, 4);

  // wait for the next loop
  t2 = millis();
  delay(fs - (t2 - t1));
}
