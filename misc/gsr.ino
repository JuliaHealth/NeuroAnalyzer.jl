const int gsr_sensor = A0;
double gsr_sensor_value = 0;
double gsr_average = 0;
double impedance = 0;

int calibration = 270;
int n = 10;
long sum = 0;

void setup(){
  Serial.begin(9600);
}

void loop(){
  sum = 0;
  // average the 10 measurements to remove the glitch
  for(int i = 0; i < n; i++) {
    gsr_sensor_value = analogRead(gsr_sensor);
    sum += gsr_sensor_value;
    delay(5);
  }
  gsr_average = sum / n;

  impedance = 1000000 * (1 / (((1024 + 2 * gsr_average) * 10000) / (calibration - gsr_average))); // convert to uSiemens
  if (impedance < 0) {
    impedance = 0;
  }
  if (isnan(impedance)) {
    impedance = 0;
  }

  Serial.print("gsr:");
  Serial.println(impedance, 4);
}
