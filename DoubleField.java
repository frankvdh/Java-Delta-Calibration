
import javax.swing.JTextField;

/**
 *
 * @author frankv
 */
public class DoubleField extends JTextField{
  double value;
  int dp;
  
public  DoubleField(double initValue, int dp) {
    super();
    this.dp = dp;
    set(initValue, dp);
  }
 
 public  DoubleField() {
    super();
    this.dp = 2;
    set(0, dp);
  }
  
  public double get() {
    try {
      value = Double.parseDouble(getText());
    } catch (NumberFormatException ex) {
      setText(fixed(value, dp));
    }
    return value;
  }
  
  public final void set(double v, int dp) {
    value = v;
    this.setText(fixed(value, dp));
  }
  
  private String fixed(double n, int dp) {
    return String.format("%." + dp + 'f', n);
  }
}
